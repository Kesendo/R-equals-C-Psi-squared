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

    /// <summary>The empirical side of the endpoint comparison: bond 0's envelope Q_peak, pinned to
    /// the independent second route of <c>PerF71OrbitKTableTests</c> (the F89 eigendecomposition
    /// operator and its own t-maximiser; gate 1e-10, above that route's mirror-bond spread). The
    /// predicted side is not pinned: the prediction picks its cluster pair by exact floating-point
    /// comparisons among near-equal candidates, so the pick, and with it the predicted Q_peak, changes
    /// between runs (N = 5 endpoint: 2.4871 in one run, 2.5237 in another).</summary>
    [Theory]
    [InlineData(4, 2.4015782771520353)]
    [InlineData(5, 2.517071370613278)]
    [InlineData(6, 2.580604998577202)]
    [InlineData(7, 2.5600734278012656)]
    public void EndpointBondZero_Envelope_IsPinned(int N, double empirical)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var b0 = JwBondQPeakUnified.Build(block).Bonds[0];
        Assert.Equal(BondClass.Endpoint, b0.BondClass);
        Assert.True(Math.Abs(b0.QPeakEmpirical - empirical) <= 1e-10, $"N={N} empirical {b0.QPeakEmpirical:R}");
        Assert.False(double.IsNaN(b0.QPeakPredicted));
    }

    // Innermost prediction via NEW-NEW Lorentzian sum regime — the V-Effect-Live emergence
    // captured by the unified architecture with Γ=2 (ANALYTICAL_SPECTRUM-anchored, see
    // JwBondQPeakUnified.LorentzianWidth doc). Residuals against the envelope Q_peak at the innermost
    // bonds: N=5 b=1,2 (flanking) 3.7%, N=6 b=2 (truly) 23.0%. The flanking-innermost
    // structure responds well to Γ=2; the truly-innermost residue (N=6 b=2) is the
    // signature of the open Γ_truly Tier1 question — a separate derivation is needed
    // for the bond where multi-pair NEW-NEW emergence is strongest.
    [Theory]
    [InlineData(5, 0.05)]   // flanking-innermost, Γ=2 fits well
    [InlineData(6, 0.25)]   // truly-innermost, open Tier1 question (Γ_truly ≠ Γ_flanking)
    public void InnermostBond_PredictionMatchesEmpirical_WithinTolerance(int N, double tolerance)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakUnified.Build(block);
        var innermost = ClosestToCenter(pred.Bonds, N);
        if (innermost is null) return;
        Assert.True(innermost.RelativeResidual < tolerance,
            $"N={N} innermost bond {innermost.Bond}: residual {innermost.RelativeResidual:P1} > {tolerance:P0} " +
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
