using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness of F50's Zeno end (Theorem E of PROOF_WEIGHT1_DEGENERACY). The exact checks hold on the
/// chain, the ring and the star, at Δ = 0, 0.1, 0.3, ½, 0.7, 1, 1.5 and 2, and each
/// fails on a wrong input fed through the same door;
/// the reading puts the slowest rate of the Liouvillian in the block the detuning form names, and its split's law is
/// the next even order.</summary>
public class ZenoEndWitnessTests
{
    private static readonly List<(int A, int B)> Chain5 = new() { (0, 1), (1, 2), (2, 3), (3, 4) };

    [Fact]
    public void Children_AreAllLive()
    {
        var w = new ZenoEndWitness(4, "chain");
        Assert.All(w.Children, c => Assert.Equal(NodeProvenance.Live, Assert.IsType<InspectableNode>(c).Provenance));
    }

    [Theory]
    [InlineData(3, "chain")]
    [InlineData(4, "ring")]
    [InlineData(8, "chain")]
    [InlineData(5, "complete")]
    public void Constructor_RejectsBadArgs(int n, string topology)
    {
        Assert.ThrowsAny<ArgumentException>(() => new ZenoEndWitness(n, topology));
    }

    [Theory]
    [InlineData(4, "chain", 1.0)]
    [InlineData(5, "chain", 1.0)]
    [InlineData(5, "chain", 0.7)]
    [InlineData(6, "chain", 0.5)]
    [InlineData(5, "chain", 0.0)]
    [InlineData(5, "ring", 1.0)]
    [InlineData(6, "ring", 2.0)]
    [InlineData(6, "ring", 0.3)]
    [InlineData(4, "star", 1.0)]
    [InlineData(5, "star", 1.5)]
    [InlineData(6, "star", 0.1)]
    public void TheExactChecks_Hold(int n, string topology, double delta)
    {
        var w = new ZenoEndWitness(n, topology, delta);
        Assert.True(w.FerromagnetHolds);
        Assert.True(w.ZenoGeneratorHolds);
        Assert.True(w.LiftHolds);
        Assert.True(w.TokenGapHolds);
        Assert.True(w.DetuningClosedFormHolds);
    }

    [Theory]
    [InlineData(5, "chain", new[] { 2, 3 })]
    [InlineData(6, "chain", new[] { 3 })]
    [InlineData(6, "ring", new[] { 3 })]
    [InlineData(5, "star", new[] { 1, 4 })]
    [InlineData(6, "star", new[] { 1, 5 })]
    public void TheFillingTheZzTermPicks_HoldsTheSlowestRate(int n, string topology, int[] expected)
    {
        var w = new ZenoEndWitness(n, topology, 1.0);
        Assert.Equal(expected, w.PredictedSlowestBlocks);
        Assert.Equal(expected, w.MeasuredSlowestBlocks);
        // the split's deviation from the detuning form falls as x^2: ratio 4 per doubling, decided at the half-integer powers
        Assert.InRange(w.DeviationRatio, Math.Pow(2, 1.5), Math.Pow(2, 2.5));
        // every block's leading coefficient tends to 2 lambda_1, its deviation falling as x^2 (ratio 4 per doubling)
        foreach (var b in w.Blocks)
            Assert.InRange(b.LeadingDeviationRatio(w.Lambda1), Math.Pow(2, 1.5), Math.Pow(2, 2.5));
    }

    [Fact]
    public void TheXyChain_TiesEveryFilling()
    {
        var w = new ZenoEndWitness(6, "chain", 0.0);
        Assert.True(w.PredictedSplitVanishes);
        Assert.True(w.ReadSplitAtRounding);
        Assert.True(w.DeviationAtRounding);
        Assert.Equal(new[] { 1, 2, 3, 4, 5 }, w.PredictedSlowestBlocks);
        Assert.Equal(new[] { 1, 2, 3, 4, 5 }, w.MeasuredSlowestBlocks);
    }

    [Fact]
    public void ASmallDelta_IsReadAsASplit_NotAsATie()
    {
        // at Δ = 0.05 the split sits far above the rounding while its deviation from δ_p − δ₁ need not: the split is
        // read, the blocks agree, and nothing says none was read
        var w = new ZenoEndWitness(5, "chain", 0.05);
        Assert.False(w.PredictedSplitVanishes);
        Assert.False(w.ReadSplitAtRounding);
        Assert.Equal(new[] { 2, 3 }, w.PredictedSlowestBlocks);
        Assert.Equal(new[] { 2, 3 }, w.MeasuredSlowestBlocks);
        Assert.DoesNotContain("none read", w.Summary, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain(w.Children, c => c.Summary.Contains("none read", StringComparison.OrdinalIgnoreCase));
    }

    [Fact]
    public void TheXyStar_TiesAtThisOrder_AndTheNextOrderPicksTheBoundary()
    {
        var w = new ZenoEndWitness(6, "star", 0.0);
        Assert.True(w.PredictedSplitVanishes);
        Assert.False(w.ReadSplitAtRounding);   // a split at x^6, above the rounding
        Assert.False(w.DeviationAtRounding);   // and nothing predicted at x^4, so the deviation is the split
        Assert.Equal(new[] { 1, 2, 3, 4, 5 }, w.PredictedSlowestBlocks);
        Assert.Equal(new[] { 1, 5 }, w.MeasuredSlowestBlocks);
        Assert.True(w.SelectionAgrees);    // the measured blocks lie in the tied set
        Assert.InRange(w.DeviationRatio, Math.Pow(2, 1.5), Math.Pow(2, 2.5));   // x^6 over x^4: x^2
    }

    [Fact]
    public void Control_ALiftAgainstARaisedSiteLaplacianFails()
    {
        var L = ZenoEndWitness.ExclusionLaplacian(5, 2, Chain5);
        var site = ZenoEndWitness.SiteLaplacian(5, Chain5);
        Assert.True(ZenoEndWitness.LiftIntertwines(L, ZenoEndWitness.LiftMatrix(5, 2), site));
        site[2, 2] += 1;
        Assert.False(ZenoEndWitness.LiftIntertwines(L, ZenoEndWitness.LiftMatrix(5, 2), site));
    }

    [Fact]
    public void Control_TheCountMovesWithADoubledBond()
    {
        var w = new ZenoEndWitness(5, "chain");
        var doubled = new List<(int A, int B)>(Chain5) { (2, 3) };
        for (int p = 1; p < 5; p++)
        {
            var L = ZenoEndWitness.ExclusionLaplacian(5, p, doubled);
            Assert.Equal(1, ZenoEndWitness.CountBelow(L, w.QHigh));   // lambda_1 moved above the old bracket
        }
    }

    [Fact]
    public void Control_TheFerromagnetAgainstAnAnisotropicMatrixFails()
    {
        var L = ZenoEndWitness.ExclusionLaplacian(5, 2, Chain5);
        var cf = ZenoEndWitness.Configurations(5, 2);
        Assert.True(ZenoEndWitness.FerromagnetIdentity(L, ZenoEndWitness.Hamiltonian(5, Chain5, 1.0).ToMatrix(), cf, 4));
        Assert.False(ZenoEndWitness.FerromagnetIdentity(L, ZenoEndWitness.Hamiltonian(5, Chain5, 2.0).ToMatrix(), cf, 4));
    }

    [Fact]
    public void Control_TheZenoGeneratorAgainstANextNearestHopFails()
    {
        var L = ZenoEndWitness.ExclusionLaplacian(5, 2, Chain5);
        var cf = ZenoEndWitness.Configurations(5, 2);
        var withHop = new List<(int A, int B)>(Chain5) { (0, 2) };
        Assert.False(ZenoEndWitness.ZenoGeneratorIsEightLaplacian(L, ZenoEndWitness.Hamiltonian(5, withHop, 1.0).ToMatrix(), cf));
    }

    [Fact]
    public void Control_TheClosedFormAgainstWeightsAtAnotherDeltaFails()
    {
        var H = ZenoEndWitness.Hamiltonian(5, Chain5, 1.0).ToMatrix();
        for (int p = 1; p < 5; p++)
        {
            var W = ZenoEndWitness.DetuningWeights(5, p, Chain5, H, ZenoEndWitness.ExactRational(0.7));
            Assert.True(ZenoEndWitness.WeightsMatchClosedForm("chain", 5, p, Chain5, W, ZenoEndWitness.ExactRational(0.7)));
            Assert.False(ZenoEndWitness.WeightsMatchClosedForm("chain", 5, p, Chain5, W, ZenoEndWitness.ExactRational(0.75)));
        }
    }

    [Fact]
    public void Control_APredictionOffByOnePercentLeavesTheWindow()
    {
        // the witness's own reading at N = 5, fed through the same deviation with every predicted split scaled by 1.01
        var w = new ZenoEndWitness(5, "chain", 1.0);
        var slowest = new double[ZenoEndWitness.ReadingXs.Length, w.Blocks.Count];
        for (int k = 0; k < ZenoEndWitness.ReadingXs.Length; k++)
            for (int p = 0; p < w.Blocks.Count; p++) slowest[k, p] = w.Blocks[p].SlowestRates[k];
        var right = w.Blocks.Select(b => b.DetuningQuotient).ToArray();
        var wrong = right.Select(q => (q - right[0]) * 1.01 + right[0]).ToArray();
        var devRight = ZenoEndWitness.SplitDeviations(slowest, right, ZenoEndWitness.ReadingXs);
        var devWrong = ZenoEndWitness.SplitDeviations(slowest, wrong, ZenoEndWitness.ReadingXs);
        Assert.InRange(devRight[1] / devRight[0], Math.Pow(2, 1.5), Math.Pow(2, 2.5));
        Assert.False(devWrong[1] / devWrong[0] > Math.Pow(2, 1.5) && devWrong[1] / devWrong[0] < Math.Pow(2, 2.5));
    }

    [Fact]
    public void Control_TheRingsFormAgainstTheChainsWeightsFails()
    {
        // a wrong count through the same door: the ring's form, every bond a bulk bond, against the chain's weights
        var H = ZenoEndWitness.Hamiltonian(5, Chain5, 1.0).ToMatrix();
        for (int p = 1; p < 5; p++)
        {
            var W = ZenoEndWitness.DetuningWeights(5, p, Chain5, H, ZenoEndWitness.ExactRational(1.0));
            Assert.True(ZenoEndWitness.WeightsMatchClosedForm("chain", 5, p, Chain5, W, ZenoEndWitness.ExactRational(1.0)));
            Assert.False(ZenoEndWitness.WeightsMatchClosedForm("ring", 5, p, Chain5, W, ZenoEndWitness.ExactRational(1.0)));
        }
    }

    [Fact]
    public void Control_EveryBulkHopTakenAsDetunedMisses()
    {
        var H = ZenoEndWitness.Hamiltonian(5, Chain5, 1.0).ToMatrix();
        var W = ZenoEndWitness.DetuningWeights(5, 3, Chain5, H, ZenoEndWitness.ExactRational(1.0));
        Assert.Equal(ZenoEndWitness.ClosedFormWeight("chain", 5, 3, 1, (1, 2)), W[(1, 2)]);
        Assert.NotEqual(ZenoEndWitness.ExactRational(2.0 * 3), W[(1, 2)]);   // 2 Delta^2 C(N-2, p-1) = 6 against 4
    }
}
