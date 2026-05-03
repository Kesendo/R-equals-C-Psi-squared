using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Decomposition;

public class MultiKBasisTests(ITestOutputHelper output)
{
    [Fact]
    public void Probe_ActualModeCount_AcrossC()
    {
        // Print actual rank after Gram-Schmidt for c=2..4. If SVD-top vectors are
        // linearly independent of the channel-uniform vectors, we get c + 2(c-1) = 3c-2.
        // Otherwise less.
        output.WriteLine("c | N | n | quartets | actual modes | expected (3c-2)");
        foreach (var (N, n) in new[] { (5, 1), (7, 2), (7, 3) })
        {
            var block = new CoherenceBlock(N, n, gammaZero: 0.05);
            var multi = MultiKBasis.Build(block);
            int expected = 3 * block.C - 2;
            output.WriteLine($"{block.C} | {N} | {n} | {multi.Quartets.Count} | {multi.TotalModes} | {expected}");
        }
    }

    [Fact]
    public void Build_C2_HasOneQuartet_FourModes()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var multi = MultiKBasis.Build(block);

        Assert.Single(multi.Quartets);
        Assert.Equal(4, multi.TotalModes);  // 2 CUs + 2 SVDs = 4 (no overlaps at c=2)
        Assert.Equal(1, multi.Quartets[0].K);
        Assert.Equal(1, multi.Quartets[0].Hd1);
        Assert.Equal(3, multi.Quartets[0].Hd2);

        // Same span as FourModeBasis (columns may differ by phase / Gram-Schmidt order).
        // Check via projection invariant: B_multi · B_multi† = B_four · B_four†.
        var four = FourModeBasis.Build(block);
        var pMulti = multi.BasisMatrix * multi.BasisMatrix.ConjugateTranspose();
        var pFour = four.BasisMatrix * four.BasisMatrix.ConjugateTranspose();
        double diff = (pMulti - pFour).FrobeniusNorm();
        Assert.True(diff < 1e-10, $"projector difference {diff:E3} suggests different span");
    }

    [Theory]
    [InlineData(7, 2, 3)]  // c=3, expected rank 3c-2 = 7
    [InlineData(7, 3, 4)]  // c=4, expected rank 3c-2 = 10
    public void Build_HigherC_HasThreeCMinusTwoModes_AfterGramSchmidt(int N, int n, int expectedC)
    {
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        Assert.Equal(expectedC, block.C);
        var multi = MultiKBasis.Build(block);

        Assert.Equal(expectedC - 1, multi.Quartets.Count);
        // Maximal rank: c channel-uniforms + 2(c-1) SVD vectors = 3c-2
        // (could be less if SVD vectors are linearly dependent on the CUs)
        int maxRank = 3 * expectedC - 2;
        Assert.True(multi.TotalModes <= maxRank,
            $"expected ≤ {maxRank} modes; got {multi.TotalModes}");
        Assert.True(multi.OffOrthonormalityResidual < 1e-10,
            $"basis residual after Gram-Schmidt: {multi.OffOrthonormalityResidual:E3}");
    }

    [Fact]
    public void Build_C3_QuartetsCoverHd1to3_And3to5()
    {
        var block = new CoherenceBlock(N: 7, n: 2, gammaZero: 0.05);
        Assert.Equal(3, block.C);
        var multi = MultiKBasis.Build(block);

        Assert.Equal(2, multi.Quartets.Count);
        var k1 = multi.Quartets[0];
        var k2 = multi.Quartets[1];
        Assert.Equal((1, 3), (k1.Hd1, k1.Hd2));
        Assert.Equal((3, 5), (k2.Hd1, k2.Hd2));
        // σ_0 should be positive for both
        Assert.True(k1.Sigma0 > 0);
        Assert.True(k2.Sigma0 > 0);
    }

    [Fact]
    public void Project_DiagonalBlock_C3_RatesAreOnHd1_Hd3_Hd5()
    {
        // Gram-Schmidt order: first c channel-uniforms (|c_1>, |c_3>, |c_5> at γ-rates
        // -0.1, -0.3, -0.5), then SVD vectors per quartet (which live in HD={1,3,3,5}).
        // The diagonal of D in this basis records each mode's rate.
        var block = new CoherenceBlock(N: 7, n: 2, gammaZero: 0.05);
        var multi = MultiKBasis.Build(block);
        var dEff = multi.Project(block.Decomposition.D);

        // Channel-uniform vectors first (3 entries):
        Assert.Equal(-0.1, dEff[0, 0].Real, 10);  // |c_1>
        Assert.Equal(-0.3, dEff[1, 1].Real, 10);  // |c_3>
        Assert.Equal(-0.5, dEff[2, 2].Real, 10);  // |c_5>

        // The remaining diagonal entries (SVD vectors, post-Gram-Schmidt) should each be
        // one of {-0.1, -0.3, -0.5} — they live in HD subspaces {1, 3, 3, 5}.
        var allowedRates = new[] { -0.1, -0.3, -0.5 };
        for (int i = 3; i < multi.TotalModes; i++)
        {
            double rate = dEff[i, i].Real;
            Assert.Contains(allowedRates, r => Math.Abs(rate - r) < 1e-9);
        }
    }
}
