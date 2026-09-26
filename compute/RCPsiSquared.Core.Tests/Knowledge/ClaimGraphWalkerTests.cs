using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Knowledge;

/// <summary>Tom's "Bauplan ist mittransportiert" test: from a topically distant
/// claim (F97, the Mandelbrot cardioid, about complex-c geometry), walking the typed
/// parent properties reaches <see cref="IF99AnchorBearing"/> foundations although F97
/// carries no F99-anchor metadata; from F99 the walk rebuilds its multi-level anchor map.
/// Alongside: genuine application-to-formula edges, the F99 diamond, and a parentless
/// formula.</summary>
public class ClaimGraphWalkerTests
{
    private static CanonicalTrigAnchorPi2Inheritance BuildF99()
    {
        var half = new HalfAsStructuralFixedPointClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var ladder = new Pi2DyadicLadderClaim();
        var staticSide = new DickeSuperpositionQuarterPi2Inheritance(ladder, quarter, half);
        var f98 = new KIntermediateAsymptoteQuarterInheritance(quarter, half, staticSide);
        return new CanonicalTrigAnchorPi2Inheritance(half, quarter, f98);
    }

    [Fact]
    public void GenuineTransitionBridgeEdge_ReachesF95()
    {
        var f95 = new F95AngleAtQuadraticZeroPi2Inheritance();
        var bridge = new TransitionBridgeF95SiblingClaim(f95);

        var reached = ClaimGraphWalker.WalkReachable(bridge);

        Assert.Equal(new Claim[] { bridge, f95 }, reached);
    }

    [Fact]
    public void GenuineLindbladSixtyDegreeEdge_ReachesF95()
    {
        var f95 = new F95AngleAtQuadraticZeroPi2Inheritance();
        var reading = new LindbladAbsorptionMatchAtSixtyDegreesClaim(f95);

        var reached = ClaimGraphWalker.WalkReachable(reading);

        Assert.Equal(new Claim[] { reading, f95 }, reached);
    }

    [Fact]
    public void F97Walk_ReachesItselfAndItsQuarterParent()
    {
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var f97 = new F97CardioidHalfFixedPointPi2Inheritance(quarter);

        var reached = ClaimGraphWalker.WalkReachable(f97);

        Assert.Equal(new Claim[] { f97, quarter }, reached);
    }

    [Fact]
    public void WalkFromF97_ReachesTheQuarterFoundation_AnF99AnchorBearingParent()
    {
        // The Bauplan test: F97 does not implement IF99AnchorBearing, but its typed
        // parent Quarter does, in the Parent role, so the walk carries the anchor
        // information to F97's vantage.
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var f97 = new F97CardioidHalfFixedPointPi2Inheritance(quarter);

        var bearing = ClaimGraphWalker.ReachableImplementing<IF99AnchorBearing>(f97);

        Assert.False(typeof(IF99AnchorBearing).IsAssignableFrom(f97.GetType()));
        var only = Assert.Single(bearing);
        Assert.Same(quarter, only);
        Assert.Equal(F99AnchorRole.Parent, only.F99Role);
    }

    [Fact]
    public void WalkFromF97_BuildsPartialMap_OneParentAndAllFiveAnchorsAsGaps()
    {
        // From F97's vantage the map is foundation-aware but holds no Direct or Covers
        // claim: nothing covered, every anchor a gap, one parent.
        var f97 = new F97CardioidHalfFixedPointPi2Inheritance(new QuarterAsBilinearMaxvalClaim());
        var claims = ClaimGraphWalker.ReachableImplementing<IF99AnchorBearing>(f97)
            .OfType<Claim>()
            .ToArray();

        var partial = new F99AnchorMap(claims);

        Assert.Empty(partial.CoveredAnchors);
        Assert.Equal(5, partial.GapAnchors.Count);
        Assert.Single(partial.ParentClaims);
    }

    [Fact]
    public void ParentlessFormula_WalksToItselfOnly()
    {
        var f95 = new F95AngleAtQuadraticZeroPi2Inheritance();

        var reached = ClaimGraphWalker.WalkReachable(f95);

        Assert.Equal(new Claim[] { f95 }, reached);
    }

    [Fact]
    public void WalkFromF99_RecursesThroughTheSixAnchorBearingClaims()
    {
        var f99 = BuildF99();

        var bearing = ClaimGraphWalker.ReachableImplementing<IF99AnchorBearing>(f99);

        Assert.Contains(bearing, b => b is CanonicalTrigAnchorPi2Inheritance);
        Assert.Contains(bearing, b => b is KIntermediateAsymptoteQuarterInheritance);
        Assert.Contains(bearing, b => b is DickeSuperpositionQuarterPi2Inheritance);
        Assert.Contains(bearing, b => b is HalfAsStructuralFixedPointClaim);
        Assert.Contains(bearing, b => b is QuarterAsBilinearMaxvalClaim);
        Assert.Contains(bearing, b => b is Pi2DyadicLadderClaim);
        Assert.Equal(6, bearing.Count);
    }

    [Fact]
    public void WalkFromF99_RebuildsItsMultiLevelAnchorMap()
    {
        var f99 = BuildF99();
        var claims = ClaimGraphWalker.ReachableImplementing<IF99AnchorBearing>(f99)
            .OfType<Claim>()
            .ToArray();

        var reconstructed = new F99AnchorMap(claims);

        Assert.Single(reconstructed.CoveredAnchors);
        Assert.Equal(3.0 / 8.0, reconstructed.CoveredAnchors[0]);
        Assert.Equal(4, reconstructed.GapAnchors.Count);
        Assert.Equal(4, reconstructed.ParentClaims.Count);
    }

    [Fact]
    public void GenuineF99Diamond_VisitsSharedParentReferencesOnce()
    {
        var f99 = BuildF99();
        Assert.Same(f99.Half, f99.F98LongTime.Half);
        Assert.Same(f99.Quarter, f99.F98LongTime.Quarter);

        var reached = ClaimGraphWalker.WalkReachable(f99);

        Assert.Single(reached, claim => ReferenceEquals(claim, f99.Half));
        Assert.Single(reached, claim => ReferenceEquals(claim, f99.Quarter));
        Assert.Equal(reached.Count, reached.Distinct(ReferenceEqualityComparer.Instance).Count());
    }

    [Fact]
    public void NullRoot_IsRejected()
    {
        Assert.Throws<ArgumentNullException>(() => ClaimGraphWalker.WalkReachable(null!));
    }
}
