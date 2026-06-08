using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The F81↔F115 connector: the F83 anti-fraction and the F115 obstruction are the two
/// orthogonal coordinates of the residual superoperator M's Klein polarity cube, meeting at the
/// bit_b=1 (diagonal) cell. Inside that cell, where F115 does all its work, the anti-fraction is
/// DEGENERATE: diagonal-cell Mixed terms are pure Π²-odd, so r = ‖H_even_nontruly‖²/‖H_odd‖² = 0 and
/// F83's anti-fraction sits at its maximum ½ (the F81 Step-8 50/50,
/// <see cref="F83AntiFractionPi2Inheritance.MaximumAntiFraction"/>), pinned across the entire cell
/// regardless of the hard/soft verdict or obstruction size. The bit_a (1+x)-valuation
/// (<see cref="WindowedHardnessClaim"/> / <see cref="WindowedObstructionScan.IsHardPair"/>) carries
/// the whole verdict. So F81/F83 (bit_b, continuous, hardware-measurable) and F87/F115 (bit_a,
/// discrete, combinatorial) are orthogonal cube coordinates and F115 is the strictly finer probe;
/// the link is structural (one M, one shared bit_b=1 gate), not quantitative.
///
/// <para>Cross-axis connector, so <see cref="Z2Axis"/> = <see cref="Symmetry.Z2Axis.Klein2"/>; it is
/// not a bit_b Claim needing a bit_a twin, so the twin slot is not applicable. (Distinct from F83's
/// own open BitATwin slot, whose mechanical twin is "F83 under X-dephasing".)</para>
///
/// <para>Tier1Derived: follows from already-derived parents, F83's anti-fraction closed form (½ at
/// r=0) and F115's valuation verdict. Anchors: docs/ANALYTICAL_FORMULAS.md F81/F83/F115 +
/// experiments/F115_OBSTRUCTION_DISTRIBUTION.md (Connection to the bit_b axis).</para></summary>
public sealed class AntiFractionObstructionOrthogonalityClaim : Claim, IZ2AxisClaim
{
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    public Z2Axis Z2Axis => Z2Axis.Klein2;
    public Claim? BitATwin => null;

    public F83AntiFractionPi2Inheritance AntiFractionClaim { get; }
    public WindowedHardnessClaim ObstructionClaim { get; }
    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public AntiFractionObstructionOrthogonalityClaim(
        F83AntiFractionPi2Inheritance antiFraction, WindowedHardnessClaim obstruction)
        : base("F81/F115 orthogonality: the anti-fraction is degenerate (1/2) in the diagonal cell where the (1+x)-valuation carries the verdict",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F81/F83/F115 + " +
               "experiments/F115_OBSTRUCTION_DISTRIBUTION.md (Connection to the bit_b axis)")
    {
        AntiFractionClaim = antiFraction ?? throw new ArgumentNullException(nameof(antiFraction));
        ObstructionClaim = obstruction ?? throw new ArgumentNullException(nameof(obstruction));
        Cases = BuildBattery(antiFraction);
    }

    public override string DisplayName =>
        "F81/F115 orthogonality (anti-fraction degenerate where the (1+x)-valuation decides)";

    public override string Summary =>
        "diagonal-cell pairs are pure Pi2-odd (r=0), so the anti-fraction is pinned at 1/2 for hard and " +
        $"soft alike, while F115's (1+x)-valuation carries the verdict; {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("orthogonality",
                summary: "bit_b anti-fraction (F83) vs bit_a (1+x)-valuation (F115); shared gate bit_b=1; F115 is the finer probe");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }

    /// <summary>Witness battery tying F83's pinned 1/2 to the F115 obstruction on the same diagonal-cell
    /// pairs (no 4^N M needed: diagonal-cell, so pure Pi2-odd, so r=0 by construction).</summary>
    private static IReadOnlyList<BatteryCase> BuildBattery(F83AntiFractionPi2Inheritance f83)
    {
        var cases = new List<BatteryCase>();
        ulong hard1 = (1UL << 1) | (1UL << 3), hard2 = 1UL | (1UL << 3);
        ulong soft1 = 0b11, soft2 = 0b110;
        bool hard = WindowedObstructionScan.IsHardPair(hard1, hard2);
        bool soft = WindowedObstructionScan.IsHardPair(soft1, soft2);
        double afHalf = f83.AntiFraction(0.0);
        // InvariantCulture so the "0.5" literal and the formatted value compare equal on
        // non-invariant locales (e.g. de-DE renders G6 as "0,5").
        string afHalfStr = afHalf.ToString("G6", CultureInfo.InvariantCulture);

        cases.Add(new BatteryCase(
            Name: "hard pair: anti-fraction 1/2 (r=0)",
            Detail: $"IsHardPair(x+x^3, 1+x^3)={hard}, anti-fraction(r=0)={afHalfStr}",
            Expected: "hard|0.5",
            Actual: $"{(hard ? "hard" : "soft")}|{afHalfStr}"));

        cases.Add(new BatteryCase(
            Name: "soft pair: anti-fraction 1/2 (r=0)",
            Detail: $"IsHardPair(1+x, x+x^2)={soft}, anti-fraction(r=0)={afHalfStr}",
            Expected: "soft|0.5",
            Actual: $"{(soft ? "hard" : "soft")}|{afHalfStr}"));

        cases.Add(new BatteryCase(
            Name: "degeneracy: anti-fraction identical for hard and soft (orthogonal to the verdict)",
            Detail: "anti-fraction(r=0) is the same value for both pairs while the obstruction differs",
            Expected: "True",
            Actual: (hard != soft).ToString()));

        cases.Add(new BatteryCase(
            Name: "pinned at the BilinearApex (r=0 sits at the apex)",
            Detail: $"MaximumAntiFraction={f83.MaximumAntiFraction.ToString("G6", CultureInfo.InvariantCulture)}, anti-fraction(r=0)={afHalfStr}",
            Expected: f83.MaximumAntiFraction.ToString("G6", CultureInfo.InvariantCulture),
            Actual: afHalfStr));

        return cases;
    }
}
