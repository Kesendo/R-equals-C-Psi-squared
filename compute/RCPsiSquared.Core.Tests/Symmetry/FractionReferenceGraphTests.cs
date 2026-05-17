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
        Assert.Equal(new[] { 1.0, 1.5, 2.0 }, FractionReferenceGraph.QBasisAnkers);
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
    public void Pi2ParityMirrors_ArePresentForOneEighthAndThreeEighths()
    {
        var g = new FractionReferenceGraph();
        // 1/8 ↔ 7/8 and 3/8 ↔ 5/8 should exist as Mirror direction
        var mirrors = g.References.Where(r => r.Direction == FractionReferenceDirection.Mirror).ToList();
        Assert.Contains(mirrors, r => Math.Abs(r.FromFraction - 1.0 / 8.0) < 1e-12
                                       && Math.Abs(r.ToFraction - 7.0 / 8.0) < 1e-12);
        Assert.Contains(mirrors, r => Math.Abs(r.FromFraction - 3.0 / 8.0) < 1e-12
                                       && Math.Abs(r.ToFraction - 5.0 / 8.0) < 1e-12);
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
