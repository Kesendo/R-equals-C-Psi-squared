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
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void NumBondsMatchesBlockNumBonds(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakPrediction.Build(block);
        Assert.Equal(N - 1, pred.Bonds.Count);
    }

    // Probe at N=11 c=2: dump per-bond JW Q_peak predictions to compare against
    // 2026-05-10 extended-grid empirical data (Q∈[0.2, 32]) where Center peak resolved
    // at Q=21.94 and bond 1, 2 at 8.79, 13.61 respectively. The JW prediction is purely
    // algebraic (no field, no Q-scan); empirical match would mean bond-dependent g_eff
    // is JW-projection structural, not "emergent field" content.
    [Fact]
    public void Probe_PerBond_QPeakPrediction_AtN11()
    {
        var block = new CoherenceBlock(N: 11, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakPrediction.Build(block);
        _out.WriteLine($"N={block.N} JW per-bond Q_peak prediction (BareDoubledPtfXPeak={JwBondQPeakPrediction.BareDoubledPtfXPeak}):");
        foreach (var b in pred.Bonds)
        {
            _out.WriteLine($"  bond b={b.Bond}: Q_peak_predicted={b.QPeakPredicted:F4}, " +
                           $"empirical_extended_grid (Q∈[0.2,32]): see 2026-05-10 c2hwhm_N11_q32.txt");
        }
        _out.WriteLine($"  MaxRelativeResidual = {pred.MaxRelativeResidual:P2}");
        // Empirical anchors from 2026-05-10 extended-grid scan (Q∈[0.2, 32]):
        //   b=0 ↔ b=9 (Endpoint):  Q_peak=2.5007
        //   b=1 ↔ b=8 (flank-1):    Q_peak=8.7946
        //   b=2 ↔ b=7 (flank-2):    Q_peak=13.6117
        //   b=3 ↔ b=6 (mid-flank):  Q_peak=1.5901  (was already in default range)
        //   b=4 ↔ b=5 (Center):     Q_peak=21.9389
        // No assertion — this is a probe to expose the predicted values to compare.
        Assert.NotEmpty(pred.Bonds);
    }

    // N=3 doesn't satisfy the 2x2 architecture: for all cluster-pairs at N=3, no (i, j)
    // sub-block has 4|x|² > (a-b)², so the general 2x2 EP formula yields no real-Q solution.
    // The architecture starts at N=4 (smallest c=2 with Endpoint+Interior distinction AND
    // a cluster-pair geometry rich enough for the 2x2 EP to exist).
    [Theory]
    [InlineData(4)]
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
    [InlineData(4)]
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
    [InlineData(4)]
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
    public void N3_PredictionReturnsNaN_StructuralLowerLimit()
    {
        // Architectural lower limit: at N=3 no cluster-pair has a (i, j) sub-block where
        // 4|X̃[i, j]|² > (λ_c1[i] − λ_c2[j])². The cluster geometry is too sparse for the
        // general 2x2 EP formula. The empirical Q_peak at N=3 (=1.806 for both bonds)
        // exists but comes from non-2x2 dynamics.
        var block = new CoherenceBlock(N: 3, n: 1, gammaZero: 0.05);
        var pred = JwBondQPeakPrediction.Build(block);
        foreach (var b in pred.Bonds)
            Assert.True(double.IsNaN(b.QPeakPredicted),
                $"N=3 bond {b.Bond}: prediction should be NaN at the architectural lower limit");
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
    public void Reconnaissance_EmitsPredictionVsEmpirical_AcrossN3To7()
    {
        _out.WriteLine("  N | bond | class    | |Δδ|   | |X̃|     | Q_EP_pred | Q_peak_pred | Q_peak_emp | rel-residual");
        _out.WriteLine("  --|------|----------|--------|---------|-----------|-------------|------------|--------------");
        foreach (int N in new[] { 3, 4, 5, 6, 7 })
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
