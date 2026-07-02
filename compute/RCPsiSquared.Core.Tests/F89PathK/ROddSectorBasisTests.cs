using RCPsiSquared.Core.F89PathK;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The R-ODD sector basis shapes (the sectorbraid deep-loci probe's substrate). The reflection
/// R splits the full (SE,DE) block into R-even and R-odd; the fixed singletons ((nBlock−1)/2 at odd
/// nBlock, 0 at even) are R-even, so oddDim = (fullDim − fixed)/2. The numeric split (full spectrum =
/// R-even ⊎ R-odd, AT q-independence) is gated in the Diagnostics tests, where the EVD lives
/// (<c>ROddDeepLociProbeTests</c>).</summary>
public class ROddSectorBasisTests
{
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void ROddBasis_HasFullDimRows_OddDimColumns_Orthonormal(int nBlock)
    {
        int fullDim = nBlock * nBlock * (nBlock - 1) / 2;                 // nBlock · C(nBlock, 2)
        int fixedSingletons = nBlock % 2 == 1 ? (nBlock - 1) / 2 : 0;
        int oddDim = (fullDim - fixedSingletons) / 2;

        var u = F89PathKSeDeBlock.ROddBasis(nBlock);
        Assert.Equal(fullDim, u.GetLength(0));
        Assert.Equal(oddDim, u.GetLength(1));

        for (int r = 0; r < oddDim; r++)                                  // UᵀU = I
            for (int s = 0; s < oddDim; s++)
            {
                double dot = 0;
                for (int i = 0; i < fullDim; i++) dot += u[i, r] * u[i, s];
                Assert.True(System.Math.Abs(dot - (r == s ? 1.0 : 0.0)) < 1e-12,
                    $"UᵀU[{r},{s}] = {dot} at nBlock={nBlock}");
            }
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void ROddAtInvariantSubspaceBasis_LivesInROddCoordinates(int k)
    {
        int nBlock = k + 1;
        int fullDim = nBlock * nBlock * (nBlock - 1) / 2;
        int fixedSingletons = nBlock % 2 == 1 ? (nBlock - 1) / 2 : 0;
        int oddDim = (fullDim - fixedSingletons) / 2;

        var basis = F89AtFactorReconstruction.ROddAtInvariantSubspaceBasis(k);
        Assert.Equal(oddDim, basis.GetLength(0));                         // rows = the R-odd sector dimension
        Assert.InRange(basis.GetLength(1), 0, oddDim);                    // cols = the R-odd AT degree
    }
}
