using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The live integer-mirror block builder, anchored on the trusted path-3 octic: the
/// ×2-cleared block's characteristic polynomial must carry the octic factor (with roots 2λ_k).</summary>
public class F89PathKSeDeBlockTests
{
    [Fact]
    public void Path3_TwoTimesBlock_HasSymDimensionTwelve()
    {
        var block = F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 2, nBlock: 4);
        Assert.Equal(12, block.GetLength(0));
        Assert.Equal(12, block.GetLength(1));
    }

    [Fact]
    public void Path3_TwoTimesBlock_CharpolyIsDivisibleByTheScaledOctic()
    {
        // The ×2-cleared block has eigenvalues 2λ_k, so its charpoly carries the octic factor with
        // roots 2λ_k = the trusted OcticCoefficientsAtQ2 literal scaled by 2^(8−k) (which is monic).
        var block = F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 2, nBlock: 4);
        var charpoly = GaussianMatrixCharpoly.Characteristic(block);

        var (re, im) = F89Path3OcticBlock.OcticCoefficientsAtQ2();
        var scaledOctic = new GaussianInteger[re.Length];
        for (int k = 0; k < re.Length; k++)
        {
            BigInteger pow = BigInteger.Pow(2, re.Length - 1 - k);   // 2^(8−k)
            scaledOctic[k] = new GaussianInteger(re[k] * pow, im[k] * pow);
        }

        var (_, remainder) = GaussianPolynomial.DivMod(charpoly, scaledOctic);
        Assert.Empty(remainder);
    }
}
