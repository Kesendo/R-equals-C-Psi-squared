using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Tests for FractionReferenceGraph: verify Tom's structural claim
/// that 0 is the universal back-reference root, all ankers connect back to
/// it, and multi-edges between fraction pairs surface the multi-viewpoint
/// structure (per Painter Principle).</summary>
public class FractionReferenceGraphTests
{
    private readonly ITestOutputHelper _out;

    public FractionReferenceGraphTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void StandardGraph_HasReferences()
    {
        var g = new FractionReferenceGraph();
        Assert.NotEmpty(g.References);
    }

    [Fact]
    public void F99Ankers_IsTheCanonicalFiveAnchorSet()
    {
        Assert.Equal(5, FractionReferenceGraph.F99Ankers.Count);
        Assert.Contains(0.0, FractionReferenceGraph.F99Ankers);
        Assert.Contains(1.0 / 8.0, FractionReferenceGraph.F99Ankers);
        Assert.Contains(1.0 / 4.0, FractionReferenceGraph.F99Ankers);
        Assert.Contains(3.0 / 8.0, FractionReferenceGraph.F99Ankers);
        Assert.Contains(1.0 / 2.0, FractionReferenceGraph.F99Ankers);
    }

    [Fact]
    public void QBasisAnkers_AreCanonicalThree()
    {
        // The 3 wave-breaking-scan anchors; now derived from QAnchorMap.CanonicalAnchors
        // by filtering for the wave_breaking_q_anchor_scan source.
        Assert.Equal(new[] { 1.0, 1.5, 2.0 }, FractionReferenceGraph.QBasisAnkers);
    }

    [Fact]
    public void NamedQAnchors_AreCanonicalTen_IncludingSqrt3()
    {
        // Full Q-anchor structure typed in QAnchorMap: onset edges, Balance, peak band
        // (with c-specific Q_peak + Q=√3 canonical θ=60°), Q_EP idealized, Endpoint.
        var named = FractionReferenceGraph.NamedQAnchors;
        Assert.Equal(10, named.Count);
        Assert.Contains(0.2, named);
        Assert.Contains(0.35, named);
        Assert.Contains(1.0, named);
        Assert.Contains(1.2, named);
        Assert.Contains(1.5, named);
        Assert.Contains(1.6, named);
        Assert.Contains(named, q => Math.Abs(q - Math.Sqrt(3.0)) < 1e-12);
        Assert.Contains(1.8, named);
        Assert.Contains(2.0, named);
        Assert.Contains(2.5, named);
    }

    [Fact]
    public void HalfToQuarter_HasMultipleViewpoints()
    {
        var g = new FractionReferenceGraph();
        var from = g.ReferencesFrom(0.5);
        var toQuarter = from.Where(r => Math.Abs(r.ToFraction - 0.25) < 1e-12).ToList();
        Assert.True(toQuarter.Count >= 4,
            $"1/2 → 1/4 should have at least 4 viewpoints (squaring + argmax/maxval + cardioid + dyadic ladder), found {toQuarter.Count}");
    }

    [Fact]
    public void F98EdgeFromThreeEighthsToQuarterIsPresent()
    {
        var g = new FractionReferenceGraph();
        var from = g.ReferencesFrom(3.0 / 8.0);
        Assert.Contains(from, r => Math.Abs(r.ToFraction - 1.0 / 4.0) < 1e-12
                                    && r.Operation.Contains("F98"));
    }

    [Fact]
    public void PolarityConvergenceEdge_AtZeroSelfLoop()
    {
        var g = new FractionReferenceGraph();
        var from = g.ReferencesFrom(0.0);
        // The (0,0) self-loop now carries Polarity direction (NOT Backward),
        // documenting the ±γ → α=0 folding convergence per
        // PolarityLayerOriginClaim. This is the structural alignment fix
        // from 2026-05-17 night agent investigation.
        Assert.Contains(from, r => Math.Abs(r.ToFraction - 0.0) < 1e-12
                                    && r.Direction == FractionReferenceDirection.Polarity
                                    && r.DocumentingClaim.Contains("XGlobalEigenstateMirror")
                                    && r.DocumentingClaim.Contains("PolarityLayerOrigin"));
    }

    [Fact]
    public void AllF99Ankers_ConvergeToMirrorAxis_TomStructuralClaim()
    {
        var g = new FractionReferenceGraph();
        Assert.True(g.AllAnkersConvergeToMirrorAxis(),
            "Tom's claim: every F99 anker has a backward-reference chain to α=0 (the polarity-mirror convergence point on the folded α-axis)");
    }

    [Fact]
    public void HalfBackwardClosure_ContainsZero()
    {
        var g = new FractionReferenceGraph();
        var closure = g.BackwardClosure(0.5);
        Assert.Contains(0.0, closure);
    }

    [Fact]
    public void EdgeCounts_ShowMultiViewpointPairs()
    {
        var g = new FractionReferenceGraph();
        var counts = g.EdgeCounts();
        // The (1/2, 1/4) pair should have ≥4 edges (the four-viewpoint test above)
        var halfQuarterCount = counts
            .Where(kv => Math.Abs(kv.Key.Item1 - 0.5) < 1e-6
                         && Math.Abs(kv.Key.Item2 - 0.25) < 1e-6)
            .Sum(kv => kv.Value);
        Assert.True(halfQuarterCount >= 4);
    }

    [Fact]
    public void Pi2ParityMirrors_AreCompleteForAllValenceComplements()
    {
        var g = new FractionReferenceGraph();
        // Full n/8 ↔ (8−n)/8 family: three non-trivial pairs + one self-mirror at n=4.
        // Periodic-table reading: alkali↔halogen, alkaline-earth↔chalcogen,
        // boron-group↔nitrogen-group, carbon-self.
        var mirrors = g.References.Where(r => r.Direction == FractionReferenceDirection.Mirror).ToList();
        Assert.Contains(mirrors, r => Math.Abs(r.FromFraction - 1.0 / 8.0) < 1e-12
                                       && Math.Abs(r.ToFraction - 7.0 / 8.0) < 1e-12);
        Assert.Contains(mirrors, r => Math.Abs(r.FromFraction - 1.0 / 4.0) < 1e-12
                                       && Math.Abs(r.ToFraction - 3.0 / 4.0) < 1e-12);
        Assert.Contains(mirrors, r => Math.Abs(r.FromFraction - 3.0 / 8.0) < 1e-12
                                       && Math.Abs(r.ToFraction - 5.0 / 8.0) < 1e-12);
    }

    [Fact]
    public void Pi2ParitySelfMirror_AtHalfAlpha_ParallelsGammaZeroPolarity()
    {
        // Carbon (n=4) is its own Π²-parity complement: 4/8 ↔ 4/8. This is
        // structurally parallel to PolarityMirrorMap's γ=0 Generic self-mirror,
        // and α=1/2 is the F86b image of γ=0 — same fixed point under
        // (1−γ²)/2 folding. Distinct from α=0's Polarity self-loop (which
        // documents ±γ → α=0 convergence, not Π²-parity).
        var g = new FractionReferenceGraph();
        var halfMirrorSelf = g.References.Where(r =>
            r.Direction == FractionReferenceDirection.Mirror
            && Math.Abs(r.FromFraction - 0.5) < 1e-12
            && Math.Abs(r.ToFraction - 0.5) < 1e-12).ToList();
        Assert.Single(halfMirrorSelf);
        Assert.Contains("carbon", halfMirrorSelf[0].Operation);
    }

    [Fact]
    public void Render_EmitsTheFullGraph()
    {
        var g = new FractionReferenceGraph();
        var rendered = g.Render();
        Assert.Contains("F99 ankers:", rendered);
        Assert.Contains("Q basis ankers:", rendered);
        Assert.Contains("AllAnkersConvergeToMirrorAxis = True", rendered);
        Assert.Contains("Multi-edge fraction pairs", rendered);
        _out.WriteLine(rendered);
    }

    [Fact]
    public void EmptyConstructor_AllowsCustomGraphs()
    {
        // A test custom graph with one edge
        var custom = new[]
        {
            new FractionReference(0.5, 0.25, "test", FractionReferenceDirection.Backward, "test"),
        };
        var g = new FractionReferenceGraph(custom);
        Assert.Single(g.References);
    }

    [Fact]
    public void NullReferences_Throws()
    {
        Assert.Throws<System.ArgumentNullException>(() =>
            new FractionReferenceGraph(null!));
    }
}
