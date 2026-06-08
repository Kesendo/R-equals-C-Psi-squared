using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The F81↔F115 connector: the F83 anti-fraction is degenerate (½) in the diagonal cell
/// where the (1+x)-valuation (F115) carries the verdict. Built directly (parents passed in) to test
/// the Claim's content independent of the registry.</summary>
public class AntiFractionObstructionOrthogonalityClaimTests
{
    private static AntiFractionObstructionOrthogonalityClaim Build()
    {
        var f83 = new F83AntiFractionPi2Inheritance(new Pi2DyadicLadderClaim());
        return new AntiFractionObstructionOrthogonalityClaim(f83, new WindowedHardnessClaim());
    }

    [Fact]
    public void Axis_IsKlein2_AndTwinSlot_NotApplicable()
    {
        var claim = Build();
        Assert.Equal(Z2Axis.Klein2, claim.Z2Axis);
        Assert.Null(claim.BitATwin);
        // BitATwinStatus is a default interface member (not overridden for Klein2, where the
        // default already yields NotApplicableForThisAxis); access it through IZ2AxisClaim, matching
        // the PolarityCubeMapTests Klein2 precedent.
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis, ((IZ2AxisClaim)claim).BitATwinStatus);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void Battery_AllCasesPass()
    {
        var claim = Build();
        Assert.True(claim.Cases.Count >= 4);
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"case '{c.Name}': expected {c.Expected}, got {c.Actual}");
    }

    /// <summary>The orthogonality itself: a hard and a soft diagonal-cell pair both have anti-fraction
    /// ½ (pure Π²-odd ⟹ r=0 ⟹ F83 apex), while their obstruction differs (hard ≥ 3, soft 0).</summary>
    [Fact]
    public void AntiFraction_IsHalf_ForBothHardAndSoft_WhileObstructionDiffers()
    {
        var f83 = new F83AntiFractionPi2Inheritance(new Pi2DyadicLadderClaim());
        ulong hard1 = (1UL << 1) | (1UL << 3), hard2 = 1UL | (1UL << 3);
        ulong soft1 = 0b11, soft2 = 0b110;
        Assert.True(WindowedObstructionScan.IsHardPair(hard1, hard2));
        Assert.False(WindowedObstructionScan.IsHardPair(soft1, soft2));
        Assert.Equal(0.5, f83.AntiFraction(0.0), 12);
        Assert.Equal(f83.MaximumAntiFraction, f83.AntiFraction(0.0), 12);
    }
}
