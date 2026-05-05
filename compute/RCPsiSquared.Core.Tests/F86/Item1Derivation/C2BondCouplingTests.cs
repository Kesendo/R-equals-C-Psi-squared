using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2BondCouplingTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void ProbeBlock_AllBondsAllEntries_MatchFourModeEffective(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        var effective = FourModeEffective.Build(block);

        for (int b = 0; b < block.NumBonds; b++)
        {
            var expected = effective.MhPerBondEff[b];
            for (int alpha = 0; alpha < 2; alpha++)
                for (int beta = 0; beta < 2; beta++)
                {
                    var actual = coupling.ProbeBlockEntry(b, alpha, beta);
                    var diff = (actual - expected[alpha, beta]).Magnitude;
                    Assert.True(diff < 1e-12,
                        $"V_b[{alpha},{beta}] at N={N}, b={b}: expected {expected[alpha, beta]}, got {actual}, diff={diff:E2}");
                }
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void ProbeBlock_OffDiagonalSumRule_F73Generalisation(int N)
    {
        // Σ_b V_b[α, β] = 0 for α ≠ β (M_H_total is diagonal in channel-uniform basis)
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);

        Complex sum01 = Complex.Zero;
        Complex sum10 = Complex.Zero;
        for (int b = 0; b < block.NumBonds; b++)
        {
            sum01 += coupling.ProbeBlockEntry(b, 0, 1);
            sum10 += coupling.ProbeBlockEntry(b, 1, 0);
        }
        Assert.True(sum01.Magnitude < 1e-12, $"Σ_b V_b[0,1] should vanish (F73), got {sum01}");
        Assert.True(sum10.Magnitude < 1e-12, $"Σ_b V_b[1,0] should vanish (F73), got {sum10}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void CrossBlock_AllBondsAllEntries_MatchFourModeEffective(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        var effective = FourModeEffective.Build(block);

        for (int b = 0; b < block.NumBonds; b++)
        {
            var expected = effective.MhPerBondEff[b];
            for (int alpha = 0; alpha < 2; alpha++)
                for (int j = 2; j < 4; j++)
                {
                    var actual = coupling.CrossBlockEntry(b, alpha, j);
                    var diff = (actual - expected[alpha, j]).Magnitude;
                    Assert.True(diff < 1e-12,
                        $"V_b[{alpha},{j}] at N={N}, b={b}: expected {expected[alpha, j]}, got {actual}, diff={diff:E2}");
                }
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void CrossBlockWitnesses_BondClassTagging_Correct(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        var witnesses = coupling.CrossBlockWitnesses;

        Assert.Equal(block.NumBonds, witnesses.Count);
        // Endpoint bonds: 0 and N-2
        Assert.Equal(BondClass.Endpoint, witnesses[0].BondClass);
        Assert.Equal(BondClass.Endpoint, witnesses[block.NumBonds - 1].BondClass);
        // Interior bonds: 1..N-3
        for (int b = 1; b < block.NumBonds - 1; b++)
            Assert.Equal(BondClass.Interior, witnesses[b].BondClass);
        // Bond index round-trip
        for (int b = 0; b < block.NumBonds; b++)
            Assert.Equal(b, witnesses[b].Bond);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void CrossBlockWitnesses_FrobeniusNorm_MatchesEntries(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);

        foreach (var w in coupling.CrossBlockWitnesses)
        {
            double expectedFrob = Math.Sqrt(
                w.EntryC1U0.Magnitude * w.EntryC1U0.Magnitude +
                w.EntryC1V0.Magnitude * w.EntryC1V0.Magnitude +
                w.EntryC3U0.Magnitude * w.EntryC3U0.Magnitude +
                w.EntryC3V0.Magnitude * w.EntryC3V0.Magnitude);
            Assert.True(Math.Abs(w.FrobeniusNorm - expectedFrob) < 1e-12,
                $"Frobenius mismatch at b={w.Bond}: stored {w.FrobeniusNorm}, recomputed {expectedFrob}");
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void CrossBlockWitness_FrobeniusNorm_BondClassSeparated(int N)
    {
        // The cross-block Frobenius norm differs between Endpoint and Interior bonds at c=2 —
        // this is the bond-position-dependence that drives the F86 Endpoint vs Interior shape
        // split. The DIRECTION of the inequality is itself an empirical observation: at
        // N=5..8 the chain-Interior cross-block Frobenius norm is larger than Endpoint
        // (e.g. N=5: endpoint=0.1237, interior=0.1934; N=8: endpoint=0.0443, interior=0.0492).
        //
        // The downstream HWHM_left/Q_peak ratio (Endpoint=0.7728 > Interior=0.7506,
        // EQ-022 (b1) 2026-05-02) is NOT a direct copy of this Frobenius ratio: it emerges
        // from the full 4×4 structure (V_b cross-block, D_eff, SVD-block) through Stage C/D.
        // Here we only assert the witness-level signature: Endpoint and Interior have
        // distinct cross-block Frobenius means.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        var witnesses = coupling.CrossBlockWitnesses;

        double endpointMean = witnesses.Where(w => w.BondClass == BondClass.Endpoint)
            .Select(w => w.FrobeniusNorm).Average();
        double interiorMean = witnesses.Where(w => w.BondClass == BondClass.Interior)
            .Select(w => w.FrobeniusNorm).Average();

        Assert.True(Math.Abs(endpointMean - interiorMean) > 1e-3,
            $"Expected endpoint vs interior cross-block Frobenius means to differ at c=2; " +
            $"got endpoint={endpointMean:F4}, interior={interiorMean:F4}");
    }

    [Fact]
    public void CrossBlockEntry_OutOfRangeAlpha_Throws()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        Assert.Throws<ArgumentOutOfRangeException>(() => coupling.CrossBlockEntry(0, alpha: 2, j: 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => coupling.CrossBlockEntry(0, alpha: -1, j: 2));
    }

    [Fact]
    public void CrossBlockEntry_OutOfRangeJ_Throws()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        Assert.Throws<ArgumentOutOfRangeException>(() => coupling.CrossBlockEntry(0, alpha: 0, j: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => coupling.CrossBlockEntry(0, alpha: 0, j: 4));
    }

    [Fact]
    public void Tier_IsTier2Verified()
    {
        // Class-level Tier reflects the weakest link: the cross-block (B2) inherits A3's
        // Tier 2 obstruction (σ_0 of V_inter exactly degenerate at even N). Probe-block
        // (B1) is structurally Tier1Derived in isolation; this is documented in the
        // class-level XML doc on C2BondCoupling.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        Assert.Equal(Tier.Tier2Verified, coupling.Tier);
    }

    [Fact]
    public void Tier_TracksInterChannelTier()
    {
        // The class Tier is sourced from the composed C2InterChannelAnalytical so future
        // promotion of A3 (e.g. closed-form ansatz reaching the 1e-10 overlap bar) flows
        // automatically into B2 without an extra hard-coded value to keep in sync.
        var block = new CoherenceBlock(N: 7, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        Assert.Equal(coupling.InterChannel.Tier, coupling.Tier);
    }

    [Fact]
    public void Build_ThrowsIfNotC2()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => C2BondCoupling.Build(block));
    }

    // ---- B3: SVD-block + AsMatrix + anti-Hermiticity guard ----------------------

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void SvdBlock_AllBondsAllEntries_MatchFourModeEffective(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        var effective = FourModeEffective.Build(block);

        for (int b = 0; b < block.NumBonds; b++)
        {
            var expected = effective.MhPerBondEff[b];
            for (int j = 2; j < 4; j++)
                for (int k = 2; k < 4; k++)
                {
                    var actual = coupling.SvdBlockEntry(b, j, k);
                    var diff = (actual - expected[j, k]).Magnitude;
                    Assert.True(diff < 1e-12,
                        $"V_b[{j},{k}] at N={N}, b={b}: expected {expected[j, k]}, got {actual}, diff={diff:E2}");
                }
        }
    }

    [Fact]
    public void SvdBlockEntry_OutOfRangeJ_Throws()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        Assert.Throws<ArgumentOutOfRangeException>(() => coupling.SvdBlockEntry(0, j: 1, k: 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => coupling.SvdBlockEntry(0, j: 4, k: 2));
    }

    [Fact]
    public void SvdBlockEntry_OutOfRangeK_Throws()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        Assert.Throws<ArgumentOutOfRangeException>(() => coupling.SvdBlockEntry(0, j: 2, k: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => coupling.SvdBlockEntry(0, j: 2, k: 4));
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void AsMatrix_FullVb_MatchesFourModeEffective(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        var effective = FourModeEffective.Build(block);

        for (int b = 0; b < block.NumBonds; b++)
        {
            var ours = coupling.AsMatrix(b);
            var theirs = effective.MhPerBondEff[b];
            var diff = (ours - theirs).FrobeniusNorm();
            Assert.True(diff < 1e-12, $"V_b at N={N}, b={b}: ‖ours − theirs‖_F = {diff:E2}");
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Vb_IsAntiHermitian_AcrossAllBondsAndEntries(int N)
    {
        // The load-bearing guard: V_b = -i [H_b, ·] projected to the 4-mode basis is
        // structurally anti-Hermitian. If probe-block, cross-block, or SVD-block accessors
        // diverge in sign convention, the three sub-blocks no longer combine to satisfy
        // the global anti-Hermiticity. This test catches such regressions.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);

        for (int b = 0; b < block.NumBonds; b++)
        {
            var Vb = coupling.AsMatrix(b);
            var anti = Vb + Vb.ConjugateTranspose();
            Assert.True(anti.FrobeniusNorm() < 1e-10,
                $"Anti-Hermiticity violated at N={N}, b={b}: ‖V_b + V_b†‖_F = {anti.FrobeniusNorm():E2}");
        }
    }

    [Fact]
    public void AsMatrix_OutOfRangeBond_Throws()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        Assert.Throws<ArgumentOutOfRangeException>(() => coupling.AsMatrix(-1));
        Assert.Throws<ArgumentOutOfRangeException>(() => coupling.AsMatrix(block.NumBonds));
    }
}
