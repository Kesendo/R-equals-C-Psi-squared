using System;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2KShapeTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void PeakOverT_MatchesFourModeResonanceScan_BitExact(int N)
    {
        // Bit-exact-equivalent check vs the existing FourModeResonanceScan path. Both must
        // compute K_b(Q, t) = 2·Re⟨ρ(t) | S_kernel | ∂ρ/∂J_b⟩ over the same (Q, t) grid and
        // pick the same peak. Tolerance 1e-10 (numerical Evd noise floor).
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var kshape = C2KShape.Build(block);
        var scan = new FourModeResonanceScan(block);
        var numerical = scan.ComputeKCurve();

        IReadOnlyList<double> qGrid = numerical.QGrid;
        IReadOnlyList<double> tGrid = numerical.TGrid;

        for (int b = 0; b < block.NumBonds; b++)
            for (int iQ = 0; iQ < qGrid.Count; iQ++)
            {
                var (kPeak, _) = kshape.PeakOverT(qGrid[iQ], b, tGrid);
                double kNumerical = numerical.KByBond[b, iQ];
                Assert.True(Math.Abs(kPeak - kNumerical) < 1e-10,
                    $"K_b(Q={qGrid[iQ]:F3}, b={b}) at N={N}: numerical={kNumerical:G6}, analytical={kPeak:G6}");
            }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void KAt_ReturnsRealValue_NotComplex(int N)
    {
        // K_b(Q, t) = 2·Re⟨ρ | S | ∂ρ⟩ is real by construction. Sanity check.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var kshape = C2KShape.Build(block);
        double k = kshape.KAt(Q: 1.5, bond: 0, t: 5.0);
        Assert.True(double.IsFinite(k));
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        // The Duhamel formula itself is exact in its inputs; the class-level Tier reflects
        // formula soundness, not the (Tier2) input quality of C2EffectiveSpectrum's
        // numerical eigenvalues.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var kshape = C2KShape.Build(block);
        Assert.Equal(Tier.Tier1Derived, kshape.Tier);
    }

    [Fact]
    public void Build_ThrowsIfNotC2()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => C2KShape.Build(block));
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void PeakOverDefaultT_AtModerateQ_PositiveAndFinite(int N)
    {
        // At Q ≈ Q_peak (Interior bond) |K| is finite and positive. Empirical Q_peak Interior
        // at c=2 N=5: 1.4821; N=7: 1.5831 — Q=1.5 sits within both bands.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var kshape = C2KShape.Build(block);
        double kPeak = kshape.PeakOverDefaultT(Q: 1.5, bond: 1);
        Assert.True(kPeak > 0);
        Assert.True(double.IsFinite(kPeak));
    }

    [Fact]
    public void KAt_OutOfRangeBond_Throws()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var kshape = C2KShape.Build(block);
        Assert.Throws<ArgumentOutOfRangeException>(() => kshape.KAt(Q: 1.0, bond: -1, t: 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => kshape.KAt(Q: 1.0, bond: block.NumBonds, t: 1.0));
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void PeakOverT_FindsTheMaximumAbsoluteK(int N)
    {
        // PeakOverT returns the envelope (max |K|, t*). An independent dense scan of the whole
        // window (2001 points, not the seed grid) must not beat it. Error model: golden section
        // stops at |Δt| ≈ √ε·t, where K lies below its maximum by ½|K''|Δt² ≈ ε·|K''|t²/2, a
        // relative shortfall of order 1e-15 here, and K itself carries round-off of the same order;
        // the gate's 1e-12 is that model with three decades to spare.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var kshape = C2KShape.Build(block);
        var tGrid = ResonanceScan.DefaultTGrid(block.GammaZero);

        const int bond = 0;
        const double Q = 1.5;
        var (peak, tAtPeak) = kshape.PeakOverT(Q, bond, tGrid);

        for (int i = 0; i <= 2000; i++)
        {
            double t = tGrid[0] + (tGrid[^1] - tGrid[0]) * i / 2000.0;
            double k = Math.Abs(kshape.KAt(Q, bond, t));
            Assert.True(k <= peak * (1 + 1e-12), $"t={t}: {k:R} above the envelope {peak:R}");
        }
        Assert.Equal(Math.Abs(kshape.KAt(Q, bond, tAtPeak)), peak, 14);
    }
}
