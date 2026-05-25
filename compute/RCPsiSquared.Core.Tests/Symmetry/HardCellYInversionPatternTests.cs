using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class HardCellYInversionPatternTests
{
    [Fact]
    public void Z2Axis_IsYParity() =>
        Assert.Equal(Z2Axis.YParity, new HardCellYInversionPattern().Z2Axis);

    [Fact]
    public void BitATwin_IsNull() =>
        Assert.Null(new HardCellYInversionPattern().BitATwin);

    [Fact]
    public void BitATwinStatus_IsNotApplicableForThisAxis() =>
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis,
            ((IZ2AxisClaim)new HardCellYInversionPattern()).BitATwinStatus);

    [Fact]
    public void Tier_IsTier1Candidate() =>
        Assert.Equal(Tier.Tier1Candidate, new HardCellYInversionPattern().Tier);
}
