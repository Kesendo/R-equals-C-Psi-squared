using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public sealed class HardCellPureDTemplateTests
{
    private readonly HardCellPureDTemplate _claim = new();

    [Fact]
    public void Tier_IsTier1Candidate() =>
        Assert.Equal(Tier.Tier1Candidate, _claim.Tier);

    [Fact]
    public void Z2Axis_IsYParity() =>
        Assert.Equal(Z2Axis.YParity, _claim.Z2Axis);

    [Fact]
    public void BitATwin_IsNull_ForYParityAxis() =>
        Assert.Null(_claim.BitATwin);

    [Fact]
    public void BitATwinStatus_IsNotApplicableForThisAxis() =>
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis,
            ((IZ2AxisClaim)_claim).BitATwinStatus);

    [Fact]
    public void Theorem_MentionsPureDTemplateAndDiagonalCell()
    {
        Assert.Contains("pure-D template", _claim.Theorem);
        Assert.Contains("diagonal Klein cell", _claim.Theorem);
        Assert.Contains("F87-hard", _claim.Theorem);
    }

    [Fact]
    public void AnchorFile_References_PROOF_F111() =>
        Assert.Contains("PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md", _claim.Anchor);

    [Fact]
    public void AnchorFile_References_F110_Parent_Observation() =>
        Assert.Contains("PROOF_F110_HARD_CELL_Y_INVERSION.md", _claim.Anchor);
}
