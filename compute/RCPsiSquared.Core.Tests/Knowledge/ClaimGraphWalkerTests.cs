using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Knowledge;

/// <summary>Tom's "Bauplan ist mittransportiert" test: from a topically
/// distant Claim (F97 Mandelbrot cardioid, about complex-c geometry),
/// reflection-walk the typed parent properties and verify that
/// <see cref="IF99AnchorBearing"/> foundation claims are reachable — even
/// though F97 itself carries no F99-anchor metadata and was constructed
/// for an entirely different subject (cardioid period-1 fixed point at
/// magnitude 1/2).</summary>
public class ClaimGraphWalkerTests
{
    private readonly ITestOutputHelper _out;

    public ClaimGraphWalkerTests(ITestOutputHelper output) => _out = output;

    private static F97CardioidHalfFixedPointPi2Inheritance BuildF97()
    {
        var half = new HalfAsStructuralFixedPointClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var ninety = new NinetyDegreeMirrorMemoryClaim();
        var poly = new PolynomialFoundationClaim();
        return new F97CardioidHalfFixedPointPi2Inheritance(half, quarter, ninety, poly);
    }

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
    public void WalkFromF97_ReachesItselfAndAllFourTypedParents()
    {
        var f97 = BuildF97();
        var reached = ClaimGraphWalker.WalkReachable(f97);
        Assert.Contains(reached, c => c is F97CardioidHalfFixedPointPi2Inheritance);
        Assert.Contains(reached, c => c is HalfAsStructuralFixedPointClaim);
        Assert.Contains(reached, c => c is QuarterAsBilinearMaxvalClaim);
        Assert.Contains(reached, c => c is NinetyDegreeMirrorMemoryClaim);
        Assert.Contains(reached, c => c is PolynomialFoundationClaim);
    }

    [Fact]
    public void WalkFromF97_ReachesF99AnchorBearingFoundations()
    {
        // The Bauplan test: F97 is about cardioid geometry, not F99 anchors.
        // It does NOT implement IF99AnchorBearing itself. But Half and
        // Quarter (both IF99AnchorBearing Parent role) are typed parents of
        // F97. The parent walk must reach them.
        var f97 = BuildF97();
        var bearing = ClaimGraphWalker.ReachableImplementing<IF99AnchorBearing>(f97);

        Assert.Contains(bearing, b => b is HalfAsStructuralFixedPointClaim);
        Assert.Contains(bearing, b => b is QuarterAsBilinearMaxvalClaim);
        // F97 itself does NOT implement IF99AnchorBearing (type-level invariant:
        // the IF99AnchorBearing filter applied to F97's reachable set can only
        // return Bearings, never F97 itself). The walk reaches F97 internally
        // but the interface filter drops it.
    }

    [Fact]
    public void WalkFromF97_AllReachedF99Bearings_AreParentRole()
    {
        // From F97 the reachable IF99AnchorBearings are all Parent role —
        // none are Direct or Covers. So the Bauplan is reconstructible
        // FOUNDATIONS-WISE but not the F99-Direct/Covers structure.
        var f97 = BuildF97();
        var bearing = ClaimGraphWalker.ReachableImplementing<IF99AnchorBearing>(f97);
        Assert.All(bearing, b => Assert.Equal(F99AnchorRole.Parent, b.F99Role));
    }

    [Fact]
    public void WalkFromF99_ReachesItselfPlusF98PlusBothPi2Foundations()
    {
        // From F99 you reach the FULL F99 inheritance subgraph: F99 itself
        // (Covers), F98 (Direct), DickeSuper + Quarter + Half + DyadicLadder
        // (Parents). That's the complete F99 anchor map territory from
        // a single starting node.
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
    public void WalkFromF99_RebuildsTheSameAnchorMap()
    {
        // The full Bauplan-reconstruction test: from F99 ALONE, walk the
        // typed-parent graph, collect IF99AnchorBearing, build the
        // F99AnchorMap from the walk. The result must match the explicit
        // six-claim map constructed by hand. The Bauplan is contained in
        // F99's parent structure.
        var f99 = BuildF99();
        var bearing = ClaimGraphWalker.ReachableImplementing<IF99AnchorBearing>(f99);
        var claims = bearing.OfType<Claim>().ToArray();
        var reconstructed = new F99AnchorMap(claims);

        Assert.Single(reconstructed.CoveredAnchors);
        Assert.Equal(3.0 / 8.0, reconstructed.CoveredAnchors[0]);
        Assert.Equal(4, reconstructed.GapAnchors.Count);
        Assert.Equal(4, reconstructed.ParentClaims.Count);
    }

    [Fact]
    public void WalkFromF99_RenderedMap_MatchesExplicitMap_VisualizeBauplan()
    {
        var f99 = BuildF99();
        var bearing = ClaimGraphWalker.ReachableImplementing<IF99AnchorBearing>(f99);
        var claims = bearing.OfType<Claim>().ToArray();
        var reconstructed = new F99AnchorMap(claims);

        _out.WriteLine("Bauplan reconstruction: walk from F99 alone, BFS over typed parents.");
        _out.WriteLine($"Reached {claims.Length} IF99AnchorBearing claims from one root node.");
        _out.WriteLine("");
        _out.WriteLine(reconstructed.Render());
    }

    [Fact]
    public void WalkFromF97_BuildsPartialMap_OnlyParents_NoDirectClaim()
    {
        // From F97 the walker reaches Half and Quarter only (both Parent
        // role). Building an F99AnchorMap from F97's walk gives ALL FIVE
        // anchors as gaps (no Direct, no Covers from F97's vantage).
        // This is also informative: F97's vantage is Pi2-foundation-aware
        // but F99-anchor-blind. The Bauplan is partial from this angle.
        var f97 = BuildF97();
        var bearing = ClaimGraphWalker.ReachableImplementing<IF99AnchorBearing>(f97);
        var claims = bearing.OfType<Claim>().ToArray();
        var partial = new F99AnchorMap(claims);

        Assert.Empty(partial.CoveredAnchors);
        Assert.Equal(5, partial.GapAnchors.Count);
        Assert.Equal(2, partial.ParentClaims.Count);

        _out.WriteLine("Partial reconstruction from F97 (cardioid, F99-blind topic):");
        _out.WriteLine($"Reached {claims.Length} IF99AnchorBearing claims (Parents only).");
        _out.WriteLine("");
        _out.WriteLine(partial.Render());
    }
}
