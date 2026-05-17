using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Tests for the X⊗N-eigenstate Mirror anchor claim: α=0 at γ=1.
/// Verifies the closed-form identity, both anchor-bearing interface
/// implementations, and that constructing the F99AnchorMap / DickeAnchorMap
/// with this claim added closes the corresponding gap.</summary>
public class XGlobalEigenstateMirrorPi2InheritanceTests
{
    private readonly ITestOutputHelper _out;

    public XGlobalEigenstateMirrorPi2InheritanceTests(ITestOutputHelper output) => _out = output;

    private static XGlobalEigenstateMirrorPi2Inheritance Build() =>
        new(new HalfAsStructuralFixedPointClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void AlphaAtMirror_IsExactlyZero()
    {
        Assert.Equal(0.0, XGlobalEigenstateMirrorPi2Inheritance.AlphaAtMirror);
    }

    [Fact]
    public void GammaAtMirror_IsExactlyOne()
    {
        Assert.Equal(1.0, XGlobalEigenstateMirrorPi2Inheritance.GammaAtMirror);
    }

    [Theory]
    [InlineData(1.0, 0.0)]    // canonical Mirror γ=+1
    [InlineData(-1.0, 0.0)]   // negative X⊗N eigenvalue γ=−1 also gives α=0
    public void AlphaFromGammaAtMirror_MatchesUniversalShape(double gamma, double expectedAlpha)
    {
        double alpha = XGlobalEigenstateMirrorPi2Inheritance.AlphaFromGammaAtMirror(gamma);
        Assert.Equal(expectedAlpha, alpha, precision: 14);
    }

    [Fact]
    public void F99Role_IsDirect()
    {
        Assert.Equal(F99AnchorRole.Direct, Build().F99Role);
    }

    [Fact]
    public void F99AnchorValues_IsSingletonZero()
    {
        var values = Build().F99AnchorValues;
        Assert.Single(values);
        Assert.Equal(0.0, values[0]);
    }

    [Fact]
    public void DickeRole_IsDirect()
    {
        Assert.Equal(DickeAnchorRole.Direct, Build().DickeRole);
    }

    [Fact]
    public void DickeAnchors_IsSingletonMirror()
    {
        var anchors = Build().DickeAnchors;
        Assert.Single(anchors);
        Assert.Equal(DickeAnchor.Mirror, anchors[0]);
    }

    [Fact]
    public void Half_TypedParent_IsExposed()
    {
        Assert.NotNull(Build().Half);
    }

    [Fact]
    public void NullHalfParent_Throws()
    {
        Assert.Throws<System.ArgumentNullException>(() =>
            new XGlobalEigenstateMirrorPi2Inheritance(null!));
    }

    [Fact]
    public void Anchor_References_DickeAnchor_And_F99_And_F98()
    {
        var anchor = Build().Anchor;
        Assert.Contains("DickeAnchor.cs", anchor);
        Assert.Contains("CanonicalTrigAnchorPi2Inheritance.cs", anchor);
        Assert.Contains("KIntermediateAsymptoteQuarterInheritance.cs", anchor);
    }

    /// <summary>Headline test: adding this claim to the F99 map closes the
    /// α=0 gap. F99 now has TWO Direct anchors (0 from Mirror, 3/8 from F98)
    /// of five, down from one. Gap count drops from 4 to 3.</summary>
    [Fact]
    public void F99AnchorMap_WithThisClaim_ClosesMirrorGap()
    {
        var half = new HalfAsStructuralFixedPointClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var ladder = new Pi2DyadicLadderClaim();
        var dickeStatic = new DickeSuperpositionQuarterPi2Inheritance(ladder, quarter, half);
        var f98 = new KIntermediateAsymptoteQuarterInheritance(quarter, half, dickeStatic);
        var f99 = new CanonicalTrigAnchorPi2Inheritance(half, quarter, f98);
        var mirror = new XGlobalEigenstateMirrorPi2Inheritance(half);

        var mapWithoutMirror = new F99AnchorMap(half, quarter, ladder, dickeStatic, f98, f99);
        var mapWithMirror = new F99AnchorMap(half, quarter, ladder, dickeStatic, f98, f99, mirror);

        // Before: 4 gaps (0, 1/8, 1/4, 1/2)
        Assert.Equal(4, mapWithoutMirror.GapAnchors.Count);
        Assert.Contains(0.0, mapWithoutMirror.GapAnchors);

        // After: 3 gaps (1/8, 1/4, 1/2) — the α=0 Mirror gap is closed
        Assert.Equal(3, mapWithMirror.GapAnchors.Count);
        Assert.DoesNotContain(0.0, mapWithMirror.GapAnchors);
        Assert.Contains(1.0 / 8.0, mapWithMirror.GapAnchors);
        Assert.Contains(1.0 / 4.0, mapWithMirror.GapAnchors);
        Assert.Contains(1.0 / 2.0, mapWithMirror.GapAnchors);

        // The new Direct claim at α=0 is the Mirror claim
        var directAtZero = mapWithMirror.DirectClaimsAt(0.0);
        Assert.Single(directAtZero);
        Assert.IsType<XGlobalEigenstateMirrorPi2Inheritance>(directAtZero[0]);

        _out.WriteLine("F99 map with Mirror Direct claim:");
        _out.WriteLine(mapWithMirror.Render());
    }

    /// <summary>Same closure on the DickeAnchor 3-set: Mirror gap closed,
    /// only Generic remains as a gap.</summary>
    [Fact]
    public void DickeAnchorMap_WithThisClaim_ClosesMirrorGap()
    {
        var half = new HalfAsStructuralFixedPointClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var ladder = new Pi2DyadicLadderClaim();
        var dickeStatic = new DickeSuperpositionQuarterPi2Inheritance(ladder, quarter, half);
        var f98 = new KIntermediateAsymptoteQuarterInheritance(quarter, half, dickeStatic);
        var f99 = new CanonicalTrigAnchorPi2Inheritance(half, quarter, f98);
        var mirror = new XGlobalEigenstateMirrorPi2Inheritance(half);

        var mapWithMirror = new DickeAnchorMap(dickeStatic, f98, f99, mirror);

        // Before: 2 gaps (Mirror, Generic); after: 1 gap (Generic only)
        Assert.Single(mapWithMirror.GapAnchors);
        Assert.Equal(DickeAnchor.Generic, mapWithMirror.GapAnchors[0]);

        var directAtMirror = mapWithMirror.DirectClaimsAt(DickeAnchor.Mirror);
        Assert.Single(directAtMirror);
        Assert.IsType<XGlobalEigenstateMirrorPi2Inheritance>(directAtMirror[0]);

        _out.WriteLine("Dicke map with Mirror Direct claim:");
        _out.WriteLine(mapWithMirror.Render());
    }
}
