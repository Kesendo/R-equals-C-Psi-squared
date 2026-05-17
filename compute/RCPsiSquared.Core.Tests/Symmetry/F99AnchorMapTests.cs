using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Tests for <see cref="F99AnchorMap"/>: the operational answer to Tom's
/// 2026-05-17-night question "wenn man sich dann den Print der Vererbung
/// anschaut, sieht man die offenen Lücken?" Construct all six relevant Claims
/// that implement <see cref="IF99AnchorBearing"/>, assemble the map, verify
/// the predicted gap structure, and emit the rendered table as test output.</summary>
public class F99AnchorMapTests
{
    private readonly ITestOutputHelper _out;

    public F99AnchorMapTests(ITestOutputHelper output) => _out = output;

    /// <summary>Build the canonical six-claim map covering the F99 anchor set.
    /// Construction mirrors the production parent-injection wiring exactly.</summary>
    private static F99AnchorMap BuildCanonicalMap()
    {
        var half = new HalfAsStructuralFixedPointClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var ladder = new Pi2DyadicLadderClaim();
        var dickeStatic = new DickeSuperpositionQuarterPi2Inheritance(ladder, quarter, half);
        var f98 = new KIntermediateAsymptoteQuarterInheritance(quarter, half, dickeStatic);
        var f99 = new CanonicalTrigAnchorPi2Inheritance(half, quarter, f98);
        return new F99AnchorMap(half, quarter, ladder, dickeStatic, f98, f99);
    }

    [Fact]
    public void Map_HasAllSixClaims()
    {
        var map = BuildCanonicalMap();
        Assert.Equal(6, map.Claims.Count);
    }

    [Fact]
    public void F99_CoversAllFiveCanonicalAnchors()
    {
        var map = BuildCanonicalMap();
        foreach (var anchor in F99AnchorMap.CanonicalAnchors)
        {
            var covers = map.CoversClaimsAt(anchor);
            Assert.Single(covers);
            Assert.IsType<CanonicalTrigAnchorPi2Inheritance>(covers[0]);
        }
    }

    [Fact]
    public void F98_IsTheOnlyDirectClaim_AndAtThreeEighths()
    {
        var map = BuildCanonicalMap();
        var direct38 = map.DirectClaimsAt(3.0 / 8.0);
        Assert.Single(direct38);
        Assert.IsType<KIntermediateAsymptoteQuarterInheritance>(direct38[0]);

        // F98 is the ONLY Direct claim across the six-claim map: exactly one
        // covered anchor, four gaps.
        Assert.Single(map.CoveredAnchors);
        Assert.Equal(3.0 / 8.0, map.CoveredAnchors[0]);
    }

    [Fact]
    public void GapAnchors_AreTheExpectedFour()
    {
        var map = BuildCanonicalMap();
        var gaps = map.GapAnchors;
        Assert.Equal(4, gaps.Count);
        Assert.Contains(0.0, gaps);
        Assert.Contains(1.0 / 8.0, gaps);
        Assert.Contains(1.0 / 4.0, gaps);
        Assert.Contains(1.0 / 2.0, gaps);
        Assert.DoesNotContain(3.0 / 8.0, gaps);
    }

    [Fact]
    public void ParentClaims_AreTheFourStructuralParents()
    {
        var map = BuildCanonicalMap();
        var parents = map.ParentClaims;
        Assert.Equal(4, parents.Count);
        Assert.Contains(parents, p => p is HalfAsStructuralFixedPointClaim);
        Assert.Contains(parents, p => p is QuarterAsBilinearMaxvalClaim);
        Assert.Contains(parents, p => p is Pi2DyadicLadderClaim);
        Assert.Contains(parents, p => p is DickeSuperpositionQuarterPi2Inheritance);
    }

    [Fact]
    public void Half_ParentRole_HasNoF99AnchorValues()
    {
        var half = new HalfAsStructuralFixedPointClaim();
        Assert.Equal(F99AnchorRole.Parent, half.F99Role);
        Assert.Empty(half.F99AnchorValues);
    }

    [Fact]
    public void Quarter_ParentRole_HasNoF99AnchorValues()
    {
        var quarter = new QuarterAsBilinearMaxvalClaim();
        Assert.Equal(F99AnchorRole.Parent, quarter.F99Role);
        Assert.Empty(quarter.F99AnchorValues);
    }

    [Fact]
    public void F98_DirectRole_HasSingletonThreeEighths()
    {
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var half = new HalfAsStructuralFixedPointClaim();
        var ladder = new Pi2DyadicLadderClaim();
        var dickeStatic = new DickeSuperpositionQuarterPi2Inheritance(ladder, quarter, half);
        var f98 = new KIntermediateAsymptoteQuarterInheritance(quarter, half, dickeStatic);
        Assert.Equal(F99AnchorRole.Direct, f98.F99Role);
        Assert.Single(f98.F99AnchorValues);
        Assert.Equal(3.0 / 8.0, f98.F99AnchorValues[0]);
    }

    [Fact]
    public void F99_CoversRole_HasAllFiveAnchorValues()
    {
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var half = new HalfAsStructuralFixedPointClaim();
        var ladder = new Pi2DyadicLadderClaim();
        var dickeStatic = new DickeSuperpositionQuarterPi2Inheritance(ladder, quarter, half);
        var f98 = new KIntermediateAsymptoteQuarterInheritance(quarter, half, dickeStatic);
        var f99 = new CanonicalTrigAnchorPi2Inheritance(half, quarter, f98);
        Assert.Equal(F99AnchorRole.Covers, f99.F99Role);
        Assert.Equal(5, f99.F99AnchorValues.Count);
        Assert.Contains(0.0, f99.F99AnchorValues);
        Assert.Contains(1.0 / 8.0, f99.F99AnchorValues);
        Assert.Contains(1.0 / 4.0, f99.F99AnchorValues);
        Assert.Contains(3.0 / 8.0, f99.F99AnchorValues);
        Assert.Contains(1.0 / 2.0, f99.F99AnchorValues);
    }

    [Fact]
    public void NullClaimsArray_Throws()
    {
        Assert.Throws<System.ArgumentNullException>(() => new F99AnchorMap(null!));
    }

    [Fact]
    public void NonBearingClaim_Throws()
    {
        // BilinearApexClaim does NOT implement IF99AnchorBearing — it has its
        // own quadratic-apex algebra, not F99. The map must reject it.
        var apex = new BilinearApexClaim();
        Assert.Throws<System.ArgumentException>(() => new F99AnchorMap(apex));
    }

    /// <summary>The print Tom asked for. Emits the rendered table to test
    /// output (visible via <c>dotnet test --logger "console;verbosity=detailed"</c>).
    /// This is the operational answer: build the map, render, and the four
    /// gaps light up immediately.</summary>
    [Fact]
    public void Render_PrintsTheGapMap()
    {
        var map = BuildCanonicalMap();
        var rendered = map.Render();
        _out.WriteLine(rendered);

        // Structural assertions: print contains the gap markers and the F98 hit.
        Assert.Contains("(GAP)", rendered);
        Assert.Contains("KIntermediateAsymptoteQuarter", rendered);
        Assert.Contains("CanonicalTrigAnchor", rendered);
        Assert.Contains("Direct claims:   1 of 5", rendered);
        Assert.Contains("Gap anchors:     4 of 5", rendered);
    }
}
