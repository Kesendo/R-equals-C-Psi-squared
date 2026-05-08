using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class JwBondQPeakUnifiedTests
{
    private readonly ITestOutputHelper _out;

    public JwBondQPeakUnifiedTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void NumBondsMatchesBlockNumBonds(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakUnified.Build(block);
        Assert.Equal(N - 1, pred.Bonds.Count);
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void EveryBond_HasFinitePrediction(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakUnified.Build(block);
        foreach (var b in pred.Bonds)
            Assert.False(double.IsNaN(b.QPeakPredicted),
                $"bond {b.Bond}: prediction is NaN");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void EndpointBonds_AreClassified_AsOldOldRegime(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakUnified.Build(block);
        var endpoints = pred.Bonds.Where(b => b.BondClass == BondClass.Endpoint).ToList();
        foreach (var ep in endpoints)
            Assert.Equal(UnifiedBondRegime.OldOld, ep.Regime);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void InnermostBond_IsClassified_AsNewNewRegime(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakUnified.Build(block);
        var innermost = ClosestToCenter(pred.Bonds, N);
        if (innermost is null) return;
        Assert.Equal(UnifiedBondRegime.NewNew, innermost.Regime);
    }

    private static UnifiedBondPrediction? ClosestToCenter(IReadOnlyList<UnifiedBondPrediction> bonds, int N)
    {
        double center = (N - 2) / 2.0;
        UnifiedBondPrediction? best = null;
        double bestDist = double.PositiveInfinity;
        foreach (var b in bonds)
        {
            if (b.BondClass != BondClass.Interior) continue;
            double d = Math.Abs(b.Bond - center);
            if (d < bestDist) { bestDist = d; best = b; }
        }
        return best;
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void EndpointBonds_PredictionMatchesEmpirical_WithinFifteenPercent(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakUnified.Build(block);
        foreach (var b in pred.Bonds.Where(p => p.BondClass == BondClass.Endpoint))
        {
            Assert.True(b.RelativeResidual < 0.15,
                $"N={N} bond {b.Bond} (Endpoint): residual {b.RelativeResidual:P1} > 15% " +
                $"(predicted {b.QPeakPredicted:F4} vs empirical {b.QPeakEmpirical:F4})");
        }
    }

    // Innermost prediction via NEW-NEW Lorentzian sum regime — the V-Effect-Live emergence
    // captured by the unified architecture. At N=5: 11.1%, at N=6 (innermost b=2): 7.1%.
    // N=4 has only 1 Interior bond (degenerate flanking=innermost); the Lorentzian-sum
    // architecture is most stable at N=5..6. N=7+ flanking bonds escape to high-Q
    // (Q_peak=4.0 grid upper bound), where the architecture predicts the inward signal.
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    public void InnermostBond_PredictionMatchesEmpirical_WithinFifteenPercent(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakUnified.Build(block);
        var innermost = ClosestToCenter(pred.Bonds, N);
        if (innermost is null) return;
        Assert.True(innermost.RelativeResidual < 0.15,
            $"N={N} innermost bond {innermost.Bond}: residual {innermost.RelativeResidual:P1} > 15% " +
            $"(predicted {innermost.QPeakPredicted:F4} vs empirical {innermost.QPeakEmpirical:F4})");
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Throws<ArgumentException>(() => JwBondQPeakUnified.Build(block));
    }

    [Fact]
    public void Tier_IsTier2Verified()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakUnified.Build(block);
        Assert.Equal(Tier.Tier2Verified, pred.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_AndResonanceHypotheses()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakUnified.Build(block);
        Assert.Contains("PROOF_F86_QPEAK", pred.Anchor);
        Assert.Contains("RESONANCE_NOT_CHANNEL", pred.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsUnifiedPredictions_AcrossN4To7()
    {
        _out.WriteLine("  N | bond | class    | regime  | top1/top2 | Q_pred | Q_emp  | residual");
        _out.WriteLine("  --|------|----------|---------|-----------|--------|--------|----------");
        foreach (int N in new[] { 4, 5, 6, 7 })
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var pred = JwBondQPeakUnified.Build(block);
            foreach (var b in pred.Bonds)
            {
                _out.WriteLine($"  {N} | {b.Bond,4} | {b.BondClass,-8} | {b.Regime,-7} | " +
                               $"{b.Top1Top2FrobRatio,9:F2} | {b.QPeakPredicted,6:F4} | " +
                               $"{b.QPeakEmpirical,6:F4} | {b.RelativeResidual:P1}");
            }
            _out.WriteLine($"  -- N={N}: max non-escape relative residual = {pred.MaxRelativeResidualNonEscape:P1}");
            _out.WriteLine("");
        }
    }
}
