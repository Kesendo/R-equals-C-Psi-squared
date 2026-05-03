using System.Numerics;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;

namespace RCPsiSquared.Core.Tests.Decomposition;

public class FourModeEffectiveTests
{
    [Theory]
    [InlineData(5, 1)]
    [InlineData(6, 1)]
    [InlineData(7, 1)]
    [InlineData(8, 1)]
    public void Build_DEff_DiagonalRates_AtC2(int N, int n)
    {
        // D in 4-mode basis: -2γ₀ on (|c_1⟩, |u_0⟩), -6γ₀ on (|c_3⟩, |v_0⟩);
        // off-diagonals zero.
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);

        Assert.Equal(-0.1, eff.DEff[0, 0].Real, 10);  // |c_1⟩
        Assert.Equal(-0.3, eff.DEff[1, 1].Real, 10);  // |c_3⟩
        Assert.Equal(-0.1, eff.DEff[2, 2].Real, 10);  // |u_0⟩
        Assert.Equal(-0.3, eff.DEff[3, 3].Real, 10);  // |v_0⟩

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                if (i != j)
                    Assert.True(eff.DEff[i, j].Magnitude < 1e-10,
                        $"D_eff[{i},{j}] should be 0 at c=2 N={N}; got {eff.DEff[i, j].Magnitude}");
    }

    [Fact]
    public void Build_ProbeEff_LivesInChannelUniformPair()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);

        Assert.True(eff.ProbeEff[0].Magnitude > 1e-3, "probe should have support on |c_1⟩");
        Assert.True(eff.ProbeEff[1].Magnitude > 1e-3, "probe should have support on |c_3⟩");
        Assert.True(eff.ProbeEff[2].Magnitude < 1e-10,
            $"probe should be ⊥ |u_0⟩; got {eff.ProbeEff[2].Magnitude}");
        Assert.True(eff.ProbeEff[3].Magnitude < 1e-10,
            $"probe should be ⊥ |v_0⟩; got {eff.ProbeEff[3].Magnitude}");
    }

    [Theory]
    [InlineData(5, 1)]
    [InlineData(6, 1)]
    [InlineData(7, 1)]
    [InlineData(8, 1)]
    public void MhTotalEff_EqualsSumOfMhPerBondEff(int N, int n)
    {
        // Linearity of projection: B† · (Σ M_b) · B = Σ B† · M_b · B
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);

        var sum = eff.MhPerBondEff[0].Clone();
        for (int b = 1; b < eff.MhPerBondEff.Count; b++) sum = sum + eff.MhPerBondEff[b];

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                Assert.True((eff.MhTotalEff[i, j] - sum[i, j]).Magnitude < 1e-10,
                    $"M_h_total_eff[{i},{j}] - Σ M_h_per_bond_eff[{i},{j}] = {(eff.MhTotalEff[i, j] - sum[i, j]).Magnitude}");
    }

    [Fact]
    public void LEffAtQ_AssemblesAsDPlusJTimesMhTotal()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);

        double q = 1.5;
        double j = q * block.GammaZero;
        var L = eff.LEffAtQ(q);

        for (int r = 0; r < 4; r++)
            for (int c = 0; c < 4; c++)
            {
                Complex expected = eff.DEff[r, c] + (Complex)j * eff.MhTotalEff[r, c];
                Assert.True((L[r, c] - expected).Magnitude < 1e-12,
                    $"L_eff[{r},{c}] mismatch at Q={q}");
            }
    }

    [Fact]
    public void MhPerBondEff_CouplesAcross2DBlocks()
    {
        // The 4-mode basis is split into [c_1, c_3] (channel-uniform / probe block)
        // and [u_0, v_0] (SVD-top / EP block). M_H_per_bond[b] should have non-trivial
        // off-block entries — that is exactly the bond-position-dependent coupling that
        // splits Endpoint vs Interior in the F86 universal-shape statement.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);

        for (int b = 0; b < eff.MhPerBondEff.Count; b++)
        {
            var Mb = eff.MhPerBondEff[b];
            double crossBlockMag = 0;
            for (int r = 0; r < 2; r++)
                for (int c = 2; c < 4; c++)
                    crossBlockMag = Math.Max(crossBlockMag, Mb[r, c].Magnitude);
            Assert.True(crossBlockMag > 1e-3,
                $"bond {b}: M_H_per_bond_eff has zero cross-block coupling — 4-mode reduction would decouple");
        }
    }
}
