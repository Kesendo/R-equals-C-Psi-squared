using RCPsiSquared.Core.CoherenceBlocks;

namespace RCPsiSquared.Core.Tests.CoherenceBlocks;

public class BlockBasisTests
{
    [Theory]
    [InlineData(5, 1, 5, 10)]   // c=2 N=5: C(5,1)·C(5,2)
    [InlineData(5, 2, 10, 10)]  // c=3 N=5
    [InlineData(7, 3, 35, 35)]  // c=4 N=7
    [InlineData(8, 1, 8, 28)]   // c=2 N=8
    [InlineData(8, 3, 56, 70)]  // c=4 N=8
    public void Dimensions_MatchBinomialCoefficients(int N, int n, int expectedMp, int expectedMq)
    {
        var basis = new BlockBasis(N, n);
        Assert.Equal(expectedMp, basis.Mp);
        Assert.Equal(expectedMq, basis.Mq);
        Assert.Equal(expectedMp * expectedMq, basis.MTotal);
    }

    [Fact]
    public void StatesP_AllHaveCorrectPopcount()
    {
        var basis = new BlockBasis(N: 5, n: 2);
        Assert.All(basis.StatesP, s =>
            Assert.Equal(2, System.Numerics.BitOperations.PopCount((ulong)s)));
        Assert.All(basis.StatesQ, s =>
            Assert.Equal(3, System.Numerics.BitOperations.PopCount((ulong)s)));
    }

    [Fact]
    public void StatesAreSortedByIntegerValue()
    {
        var basis = new BlockBasis(N: 6, n: 2);
        for (int i = 1; i < basis.StatesP.Count; i++)
            Assert.True(basis.StatesP[i] > basis.StatesP[i - 1]);
        for (int i = 1; i < basis.StatesQ.Count; i++)
            Assert.True(basis.StatesQ[i] > basis.StatesQ[i - 1]);
    }

    [Fact]
    public void FlatIndex_IsConsistent()
    {
        var basis = new BlockBasis(N: 5, n: 1);
        // For (p=1, q=3) at N=5: p has popcount 1 (bit 0 set, big-endian → site 4),
        // q has popcount 2. FlatIndex = IndexP(1)·Mq + IndexQ(3).
        long p = basis.StatesP[2]; // some valid state
        long q = basis.StatesQ[5];
        int idx = basis.FlatIndex(p, q);
        Assert.Equal(basis.IndexP(p) * basis.Mq + basis.IndexQ(q), idx);
        Assert.InRange(idx, 0, basis.MTotal - 1);
    }

    [Fact]
    public void IndexP_AndIndexQ_AreInverseOfStatesLookup()
    {
        var basis = new BlockBasis(N: 6, n: 2);
        for (int i = 0; i < basis.Mp; i++)
            Assert.Equal(i, basis.IndexP(basis.StatesP[i]));
        for (int i = 0; i < basis.Mq; i++)
            Assert.Equal(i, basis.IndexQ(basis.StatesQ[i]));
    }
}
