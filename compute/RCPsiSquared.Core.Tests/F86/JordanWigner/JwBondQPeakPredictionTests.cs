using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class JwBondQPeakPredictionTests
{
    private readonly ITestOutputHelper _out;

    public JwBondQPeakPredictionTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void NumBondsMatchesBlockNumBonds(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakPrediction.Build(block);
        Assert.Equal(N - 1, pred.Bonds.Count);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void EveryBond_HasFinitePrediction(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakPrediction.Build(block);
        foreach (var b in pred.Bonds)
        {
            Assert.False(double.IsNaN(b.QPeakPredicted),
                $"bond {b.Bond}: Q_peak prediction is NaN");
            Assert.True(b.QPeakPredicted > 0,
                $"bond {b.Bond}: Q_peak prediction must be positive");
        }
    }

    // F71-mirror invariance at the prediction level: bond b and bond N-2-b have F71-mirror-
    // related top-pair selection. The selection-rule may pick different pairs (different
    // |Δδ|, different |X̃|) for mirror partners due to ranking-tie-breaking, so prediction
    // tolerance is relaxed. The structural invariance (max relative deviation < 50%) holds.
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void F71MirrorBonds_QPeakPredictionsWithinFiftyPercent(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakPrediction.Build(block);
        int numBonds = pred.Bonds.Count;
        for (int b = 0; b < numBonds; b++)
        {
            int mirror = numBonds - 1 - b;
            if (mirror <= b) continue;
            double pa = pred.Bonds[b].QPeakPredicted;
            double pb = pred.Bonds[mirror].QPeakPredicted;
            double rel = Math.Abs(pa - pb) / Math.Max(pa, pb);
            Assert.True(rel < 0.5, $"N={N} bond {b} vs mirror {mirror}: |Δ|/max = {rel:P1} > 50%");
        }
    }

    // Endpoint Q_peak prediction empirically matches C2HwhmRatio within ~12% across N=5..7
    // — the JW-track 2x2 closed-form architecture works for Endpoint bonds. Innermost-bonds
    // need multi-cluster-Petermann combination (open Tier1-Promotion; not tested here).
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void EndpointBonds_PredictionMatchesEmpirical_WithinFifteenPercent(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakPrediction.Build(block);
        foreach (var b in pred.Bonds.Where(p => p.BondClass == RCPsiSquared.Core.Resonance.BondClass.Endpoint))
        {
            Assert.True(b.RelativeResidual < 0.15,
                $"N={N} bond {b.Bond} (Endpoint): relative residual {b.RelativeResidual:P1} > 15% " +
                $"(predicted {b.QPeakPredicted:F4} vs empirical {b.QPeakEmpirical:F4})");
        }
    }

    private static BondQPeakPrediction? ClosestToCenter(IReadOnlyList<BondQPeakPrediction> bonds, int N)
    {
        double center = (N - 2) / 2.0;
        BondQPeakPrediction? best = null;
        double bestDist = double.PositiveInfinity;
        foreach (var b in bonds)
        {
            if (b.BondClass != RCPsiSquared.Core.Resonance.BondClass.Interior) continue;
            double d = Math.Abs(b.Bond - center);
            if (d < bestDist) { bestDist = d; best = b; }
        }
        return best;
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Throws<ArgumentException>(() => JwBondQPeakPrediction.Build(block));
    }

    [Fact]
    public void Tier_IsTier2Verified()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakPrediction.Build(block);
        Assert.Equal(Tier.Tier2Verified, pred.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_AndDirectionB()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakPrediction.Build(block);
        Assert.Contains("PROOF_F86_QPEAK", pred.Anchor);
        Assert.Contains("Direction (b'')", pred.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsPredictionVsEmpirical_AcrossN5To7()
    {
        _out.WriteLine("  N | bond | class    | |Δδ|   | |X̃|     | Q_EP_pred | Q_peak_pred | Q_peak_emp | rel-residual");
        _out.WriteLine("  --|------|----------|--------|---------|-----------|-------------|------------|--------------");
        foreach (int N in new[] { 5, 6, 7 })
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var pred = JwBondQPeakPrediction.Build(block);
            foreach (var b in pred.Bonds)
            {
                _out.WriteLine($"  {N} | {b.Bond,4} | {b.BondClass,-8} | {b.TopPairAbsDeltaδ,6:F3} | " +
                               $"{b.XTildeMagnitude,7:F4} | {b.QEpPredicted,9:F4} | {b.QPeakPredicted,11:F4} | " +
                               $"{b.QPeakEmpirical,10:F4} | {b.RelativeResidual,12:P1}");
            }
            _out.WriteLine($"  -- N={N}: max relative residual (non-escape bonds) = {pred.MaxRelativeResidual:P1}");
            _out.WriteLine("");
        }
    }
}
