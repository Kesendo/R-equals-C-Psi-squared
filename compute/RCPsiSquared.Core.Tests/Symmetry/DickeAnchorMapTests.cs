using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Tests for <see cref="DickeAnchorMap"/>: the second I*Bearing+*Map
/// pair (mirror of <see cref="F99AnchorMapTests"/>). The DickeAnchor 3-set
/// {Mirror, KIntermediate, Generic} is a structural subset of F99's 5-set;
/// this map surfaces gap structure specifically through the uniform-Dicke
/// lens.</summary>
public class DickeAnchorMapTests
{
    private readonly ITestOutputHelper _out;

    public DickeAnchorMapTests(ITestOutputHelper output) => _out = output;

    private static DickeAnchorMap BuildCanonicalMap()
    {
        var half = new HalfAsStructuralFixedPointClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var ladder = new Pi2DyadicLadderClaim();
        var dickeStatic = new DickeSuperpositionQuarterPi2Inheritance(ladder, quarter, half);
        var f98 = new KIntermediateAsymptoteQuarterInheritance(quarter, half, dickeStatic);
        var f99 = new CanonicalTrigAnchorPi2Inheritance(half, quarter, f98);
        // The three IDickeAnchorBearing claims (F98 Direct, DickeSuper Direct, F99 Covers).
        return new DickeAnchorMap(dickeStatic, f98, f99);
    }

    [Fact]
    public void Map_HasThreeBearingClaims()
    {
        var map = BuildCanonicalMap();
        Assert.Equal(3, map.Claims.Count);
    }

    [Fact]
    public void F99_CoversAllThreeDickeAnchors()
    {
        var map = BuildCanonicalMap();
        foreach (var anchor in DickeAnchorMap.CanonicalAnchors)
        {
            var covers = map.CoversClaimsAt(anchor);
            Assert.Single(covers);
            Assert.IsType<CanonicalTrigAnchorPi2Inheritance>(covers[0]);
        }
    }

    [Fact]
    public void KIntermediate_HasTwoDirectClaims_F98AndDickeStatic()
    {
        var map = BuildCanonicalMap();
        var direct = map.DirectClaimsAt(DickeAnchor.KIntermediate);
        Assert.Equal(2, direct.Count);
        Assert.Contains(direct, c => c is KIntermediateAsymptoteQuarterInheritance);
        Assert.Contains(direct, c => c is DickeSuperpositionQuarterPi2Inheritance);
    }

    [Fact]
    public void MirrorAndGeneric_AreGaps()
    {
        var map = BuildCanonicalMap();
        Assert.Empty(map.DirectClaimsAt(DickeAnchor.Mirror));
        Assert.Empty(map.DirectClaimsAt(DickeAnchor.Generic));
        Assert.Equal(2, map.GapAnchors.Count);
        Assert.Contains(DickeAnchor.Mirror, map.GapAnchors);
        Assert.Contains(DickeAnchor.Generic, map.GapAnchors);
    }

    [Fact]
    public void OnlyKIntermediate_IsCovered_OfThreeAnchors()
    {
        var map = BuildCanonicalMap();
        Assert.Single(map.CoveredAnchors);
        Assert.Equal(DickeAnchor.KIntermediate, map.CoveredAnchors[0]);
    }

    [Fact]
    public void F98_DickeRole_IsDirect_AtKIntermediate()
    {
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var half = new HalfAsStructuralFixedPointClaim();
        var ladder = new Pi2DyadicLadderClaim();
        var dickeStatic = new DickeSuperpositionQuarterPi2Inheritance(ladder, quarter, half);
        var f98 = new KIntermediateAsymptoteQuarterInheritance(quarter, half, dickeStatic);
        Assert.Equal(DickeAnchorRole.Direct, f98.DickeRole);
        Assert.Single(f98.DickeAnchors);
        Assert.Equal(DickeAnchor.KIntermediate, f98.DickeAnchors[0]);
    }

    [Fact]
    public void DickeStatic_DickeRole_IsDirect_AtKIntermediate()
    {
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var half = new HalfAsStructuralFixedPointClaim();
        var ladder = new Pi2DyadicLadderClaim();
        var dickeStatic = new DickeSuperpositionQuarterPi2Inheritance(ladder, quarter, half);
        Assert.Equal(DickeAnchorRole.Direct, dickeStatic.DickeRole);
        Assert.Single(dickeStatic.DickeAnchors);
        Assert.Equal(DickeAnchor.KIntermediate, dickeStatic.DickeAnchors[0]);
    }

    [Fact]
    public void F99_DickeRole_IsCovers_AllThreeAnchors()
    {
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var half = new HalfAsStructuralFixedPointClaim();
        var ladder = new Pi2DyadicLadderClaim();
        var dickeStatic = new DickeSuperpositionQuarterPi2Inheritance(ladder, quarter, half);
        var f98 = new KIntermediateAsymptoteQuarterInheritance(quarter, half, dickeStatic);
        var f99 = new CanonicalTrigAnchorPi2Inheritance(half, quarter, f98);
        Assert.Equal(DickeAnchorRole.Covers, f99.DickeRole);
        Assert.Equal(3, f99.DickeAnchors.Count);
        Assert.Contains(DickeAnchor.Mirror, f99.DickeAnchors);
        Assert.Contains(DickeAnchor.KIntermediate, f99.DickeAnchors);
        Assert.Contains(DickeAnchor.Generic, f99.DickeAnchors);
    }

    [Fact]
    public void NonBearingClaim_Throws()
    {
        var apex = new BilinearApexClaim();
        Assert.Throws<System.ArgumentException>(() => new DickeAnchorMap(apex));
    }

    /// <summary>The print: emits the rendered DickeAnchor map to test output.
    /// Shows Mirror gap and Generic gap explicitly; KIntermediate has two
    /// Direct claims (F98 dynamic-side + DickeSuper static-side).</summary>
    [Fact]
    public void Render_PrintsTheDickeGapMap()
    {
        var map = BuildCanonicalMap();
        var rendered = map.Render();
        _out.WriteLine(rendered);

        Assert.Contains("(GAP)", rendered);
        Assert.Contains("Mirror", rendered);
        Assert.Contains("KIntermediate", rendered);
        Assert.Contains("Generic", rendered);
        Assert.Contains("KIntermediateAsymptoteQuarter", rendered);
        Assert.Contains("DickeSuperpositionQuarter", rendered);
        Assert.Contains("Direct claims:   1 of 3", rendered);
        Assert.Contains("Gap anchors:     2 of 3", rendered);
    }
}
