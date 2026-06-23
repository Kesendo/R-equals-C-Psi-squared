using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The path-3 anchor for full Option D: isolate F_d LIVE by dividing the reconstructed AT
/// factor out of the block's characteristic polynomial (validation triple), and confirm the
/// result equals the trusted octic — proving the isolated polynomial IS the H_B-mixed factor,
/// not an imported literal.</summary>
public class F89HbMixedIsolationTests
{
    private static GaussianInteger[] ScaledOctic()
    {
        var (re, im) = F89Path3OcticBlock.OcticCoefficientsAtQ2();
        var scaled = new GaussianInteger[re.Length];
        for (int k = 0; k < re.Length; k++)
        {
            BigInteger pow = BigInteger.Pow(2, re.Length - 1 - k);   // 2^(8−k): roots 2λ_k
            scaled[k] = new GaussianInteger(re[k] * pow, im[k] * pow);
        }
        return scaled;
    }

    [Fact]
    public void Path3_IsolatesFdViaAtFactor_AndItEqualsTheScaledOctic()
    {
        var block = F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 2, nBlock: 4);
        var charpoly = GaussianMatrixCharpoly.Characteristic(block);
        var atScaled = F89AtFactorReconstruction.ForPath3(q0: 2);

        var fd = F89HbMixedIsolation.Isolate(charpoly, atScaled, expectedDegree: 8);

        Assert.Equal(ScaledOctic(), fd);
    }

    [Fact]
    public void Isolate_WithAWrongAtFactor_ThrowsInsteadOfReturningGarbage()
    {
        var block = F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 2, nBlock: 4);
        var charpoly = GaussianMatrixCharpoly.Characteristic(block);

        // (λ−1)(λ−2)(λ−3)(λ−5): a monic degree-4 polynomial that does NOT divide the charpoly.
        var wrongAt = GaussianPolynomial.Multiply(
            GaussianPolynomial.Multiply(new GaussianInteger[] { -1, 1 }, new GaussianInteger[] { -2, 1 }),
            GaussianPolynomial.Multiply(new GaussianInteger[] { -3, 1 }, new GaussianInteger[] { -5, 1 }));

        Assert.Throws<System.InvalidOperationException>(
            () => F89HbMixedIsolation.Isolate(charpoly, wrongAt, expectedDegree: 8));
    }
}
