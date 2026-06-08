using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The "two blind spots of the anti-fraction" connector, built directly (parents passed in).</summary>
public class AntiFractionTwoBlindSpotsClaimTests
{
    private static AntiFractionTwoBlindSpotsClaim Build()
    {
        var f83 = new F83AntiFractionPi2Inheritance(new Pi2DyadicLadderClaim());
        var f89 = new F89F87BreakPredictionFromF83(new F89F87TrulyInheritance(new F87TrichotomyClassification()));
        var halfEnd = new AntiFractionObstructionOrthogonalityClaim(f83, new WindowedHardnessClaim());
        return new AntiFractionTwoBlindSpotsClaim(f89, halfEnd);
    }

    [Fact]
    public void Axis_IsKlein2_AndTwinSlot_NotApplicable()
    {
        var claim = Build();
        Assert.Equal(Z2Axis.Klein2, claim.Z2Axis);
        Assert.Null(claim.BitATwin);
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis, ((IZ2AxisClaim)claim).BitATwinStatus);
    }

    [Fact]
    public void Tier_IsTier1Derived() => Assert.Equal(Tier.Tier1Derived, Build().Tier);

    [Fact]
    public void Battery_AllCasesPass()
    {
        var claim = Build();
        Assert.True(claim.Cases.Count >= 4);
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"case '{c.Name}': expected {c.Expected}, got {c.Actual}");
    }
}
