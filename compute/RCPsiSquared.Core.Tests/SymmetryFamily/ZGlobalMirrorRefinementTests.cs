using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.SymmetryFamily;
using Xunit;

namespace RCPsiSquared.Core.Tests.SymmetryFamily;

public sealed class ZGlobalMirrorRefinementTests
{
    [Fact]
    public void Tier_IsTier1Derived()
    {
        var z = new ZGlobalMirrorRefinement(new SymmetryFamilyInventory());
        Assert.Equal(Tier.Tier1Derived, z.Tier);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void ParityFromJointPopcount_MatchesZGlobalEigenvalue(int N)
    {
        for (int pCol = 0; pCol <= N; pCol++)
            for (int pRow = 0; pRow <= N; pRow++)
            {
                int parityFromJointPopcount = (pCol + pRow) % 2;
                int parityFromZGlobal = ZGlobalMirrorRefinement.JointPopcountParity(pCol, pRow);
                Assert.Equal(parityFromJointPopcount, parityFromZGlobal);
            }
    }

    [Fact]
    public void Summary_MentionsRedundancyWithJointPopcount()
    {
        var z = new ZGlobalMirrorRefinement(new SymmetryFamilyInventory());
        Assert.Contains("redundant", z.Summary, System.StringComparison.OrdinalIgnoreCase);
    }
}
