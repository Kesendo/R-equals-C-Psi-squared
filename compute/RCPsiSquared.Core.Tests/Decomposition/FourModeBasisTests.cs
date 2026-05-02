using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Probes;

namespace RCPsiSquared.Core.Tests.Decomposition;

public class FourModeBasisTests
{
    [Theory]
    [InlineData(5, 1)]
    [InlineData(6, 1)]
    [InlineData(7, 1)]
    [InlineData(8, 1)]
    public void Basis_IsOrthonormal_AtC2(int N, int n)
    {
        // step_i finding: {|c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩} are mutually orthonormal at c=2.
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var basis = FourModeBasis.Build(block);
        Assert.True(basis.OffOrthonormalityResidual < 1e-10,
            $"4-mode basis not orthonormal at c=2 N={N}: max |off-diag| = {basis.OffOrthonormalityResidual}");
    }

    [Fact]
    public void Project_DickeProbe_HasSupportOnlyInChannelUniformPair()
    {
        // probe has support on |c_1⟩ and |c_3⟩ but NOT on |u_0⟩, |v_0⟩.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var basis = FourModeBasis.Build(block);
        var probe = DickeBlockProbe.Build(block);
        var probe4 = basis.Project(probe);

        Assert.True(probe4[0].Magnitude > 1e-3,
            "probe should have support on |c_1⟩");
        Assert.True(probe4[1].Magnitude > 1e-3,
            "probe should have support on |c_3⟩");
        Assert.True(probe4[2].Magnitude < 1e-10,
            $"probe should be ⊥ |u_0⟩; got {probe4[2].Magnitude}");
        Assert.True(probe4[3].Magnitude < 1e-10,
            $"probe should be ⊥ |v_0⟩; got {probe4[3].Magnitude}");
    }

    [Fact]
    public void Project_BlockL_DiagonalIsExpectedRates()
    {
        // D in 4-mode basis: -2γ₀ on (|c_1⟩, |u_0⟩), -6γ₀ on (|c_3⟩, |v_0⟩).
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var basis = FourModeBasis.Build(block);
        var d4 = basis.Project(block.Decomposition.D);

        Assert.Equal(-0.1, d4[0, 0].Real, 10);  // |c_1⟩
        Assert.Equal(-0.3, d4[1, 1].Real, 10);  // |c_3⟩
        Assert.Equal(-0.1, d4[2, 2].Real, 10);  // |u_0⟩
        Assert.Equal(-0.3, d4[3, 3].Real, 10);  // |v_0⟩

        // D should be diagonal in this basis (no cross-channel mixing within D)
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                if (i == j) continue;
                Assert.True(d4[i, j].Magnitude < 1e-10,
                    $"D[{i},{j}] should be 0; got {d4[i, j].Magnitude}");
            }
        }
    }
}
