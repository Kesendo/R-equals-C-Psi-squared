using System.Numerics;
using RCPsiSquared.Core.CoherenceBlocks;
using Xunit;

namespace RCPsiSquared.Core.Tests.CoherenceBlocks;

public class BlockLColumnMajorRawTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void BuildUniformLColumnMajorRaw_AgreesWithBuildUniformLAt(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        double j = 0.075;
        var matrixForm = BlockLDecomposition.BuildUniformLAt(block, j);
        Complex[] raw = BlockLDecomposition.BuildUniformLColumnMajorRaw(block, j);

        int mTotal = block.Basis.MTotal;
        Assert.Equal((long)mTotal * mTotal, raw.LongLength);
        for (int row = 0; row < mTotal; row++)
            for (int col = 0; col < mTotal; col++)
            {
                Complex expected = matrixForm[row, col];
                Complex actual = raw[(long)col * mTotal + row];   // column-major: idx = col*mTotal + row
                Assert.True((expected - actual).Magnitude < 1e-12,
                    $"L[{row},{col}]: matrix={expected}, raw={actual}");
            }
    }
}
