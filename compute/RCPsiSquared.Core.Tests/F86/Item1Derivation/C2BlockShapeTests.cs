using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2BlockShapeTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void BlockDimensions_C2_MatchAnalyticalCounts(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var shape = new C2BlockShape(block);
        Assert.Equal(2, block.C);
        Assert.Equal(N, shape.PnDimension);
        Assert.Equal(N * (N - 1) / 2, shape.PnPlus1Dimension);
        Assert.Equal(N * (N - 1), shape.HdEqualsOnePairs);
        Assert.Equal(N * (N - 1) * (N - 2) / 2, shape.HdEqualsThreePairs);
        Assert.Equal(Tier.Tier1Derived, shape.Tier);
    }

    [Fact]
    public void Constructor_ThrowsIfNotC2()
    {
        // c=3 block: N=5, n=2 gives c = min(2, 2) + 1 = 3
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Equal(3, block.C);
        Assert.Throws<ArgumentException>(() => new C2BlockShape(block));
    }
}
