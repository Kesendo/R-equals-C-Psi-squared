using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Probes;

namespace RCPsiSquared.Core.Tests.Probes;

public class SpatialSumKernelTests
{
    [Theory]
    [InlineData(5, 1)]
    [InlineData(7, 3)]
    [InlineData(6, 2)]
    public void Kernel_IsHermitian(int N, int n)
    {
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var s = SpatialSumKernel.Build(block);

        for (int i = 0; i < s.RowCount; i++)
        {
            for (int j = 0; j < s.ColumnCount; j++)
            {
                Assert.Equal(s[i, j].Real, s[j, i].Real, 12);
                Assert.Equal(s[i, j].Imaginary, -s[j, i].Imaginary, 12);
            }
        }
    }

    [Theory]
    [InlineData(5, 1)]
    [InlineData(6, 2)]
    public void Kernel_IsPositiveSemiDefinite(int N, int n)
    {
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var s = SpatialSumKernel.Build(block);

        var dense = s.ToArray();
        var herm = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>.Build.DenseOfArray(dense);
        var evd = herm.Evd();
        for (int i = 0; i < evd.EigenValues.Count; i++)
        {
            double lambda = evd.EigenValues[i].Real;
            Assert.True(lambda > -1e-10, $"S has eigenvalue {lambda} < 0; not PSD.");
            Assert.True(Math.Abs(evd.EigenValues[i].Imaginary) < 1e-10,
                $"S has complex eigenvalue {evd.EigenValues[i]}; not Hermitian.");
        }
    }

    [Fact]
    public void F73_SpatialSumPurity_AtZero_IsOneHalf_For_C1_Block()
    {
        // F73: at c=1 (vac-SE / (0,1) block), the Dicke probe |vac⟩⟨S_1|/2 has S(t=0) = 1/2 exactly.
        var block = new CoherenceBlock(N: 4, n: 0, gammaZero: 0.05);
        Assert.Equal(1, block.C);
        var probe = DickeBlockProbe.Build(block);
        var s = SpatialSumKernel.Build(block);

        var sProbe = s * probe;
        var quad = probe.ConjugateDotProduct(sProbe);
        Assert.Equal(0.5, quad.Real, 12);
        Assert.Equal(0.0, quad.Imaginary, 12);
    }
}
