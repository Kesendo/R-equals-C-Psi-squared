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
}
