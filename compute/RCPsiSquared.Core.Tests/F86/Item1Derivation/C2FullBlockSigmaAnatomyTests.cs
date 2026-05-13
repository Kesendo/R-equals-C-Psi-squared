using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2FullBlockSigmaAnatomyTests
{
    private static CoherenceBlock C2Block(int N) =>
        new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);

    [Fact]
    public void Build_AtC2N5_ReturnsTier2VerifiedClaim()
    {
        var anatomy = C2FullBlockSigmaAnatomy.Build(C2Block(5));
        Assert.Equal(Tier.Tier2Verified, anatomy.Tier);
        Assert.Equal(5, anatomy.Block.N);
    }

    [Fact]
    public void Build_AtC2N5_ProducesOneWitnessPerEigenmode()
    {
        var block = C2Block(5);
        var anatomy = C2FullBlockSigmaAnatomy.Build(block);
        Assert.Equal(block.Basis.MTotal, anatomy.SigmaSpectrum.Count);
    }

    [Fact]
    public void SigmaWitness_AtC2N5_HasNonNegativeSigma()
    {
        var anatomy = C2FullBlockSigmaAnatomy.Build(C2Block(5));
        Assert.All(anatomy.SigmaSpectrum, w =>
            Assert.True(w.Sigma >= -1e-12,
                $"Sigma must be non-negative (S is PSD); got {w.Sigma} at λ={w.EigenvalueReal}+{w.EigenvalueImag}i"));
    }

    [Fact]
    public void SigmaWitness_AtC2N5_HasPositiveTotalSigma()
    {
        var anatomy = C2FullBlockSigmaAnatomy.Build(C2Block(5));
        double total = anatomy.SigmaSpectrum.Sum(w => w.Sigma);
        Assert.True(total > 1e-6, $"Total sigma must be positive; got {total}");
    }

    [Theory]
    [InlineData(5, 2)]   // c=2 N=5: path-4 → F_a count = floor(5/2) = 2
    [InlineData(6, 3)]   // c=2 N=6: path-5 → F_a count = floor(6/2) = 3
    [InlineData(7, 3)]   // c=2 N=7: path-6 → F_a count = floor(7/2) = 3
    [InlineData(8, 4)]   // c=2 N=8: path-7 → F_a count = floor(8/2) = 4
    public void FaModes_Count_MatchesFaCount(int n, int expectedFaCount)
    {
        var anatomy = C2FullBlockSigmaAnatomy.Build(C2Block(n));
        int actualFaCount = anatomy.SigmaSpectrum.Count(w => w.BlochIndexN.HasValue);
        Assert.Equal(expectedFaCount, actualFaCount);
    }


    [Fact]
    public void FaModes_BlochIndices_AreInSeAntiOrbit()
    {
        var anatomy = C2FullBlockSigmaAnatomy.Build(C2Block(7));   // N=7 → orbit {2, 4, 6}
        var assigned = anatomy.SigmaSpectrum
            .Where(w => w.BlochIndexN.HasValue)
            .Select(w => w.BlochIndexN!.Value)
            .OrderBy(n => n)
            .ToArray();
        Assert.Equal(new[] { 2, 4, 6 }, assigned);
    }
}
