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
}
