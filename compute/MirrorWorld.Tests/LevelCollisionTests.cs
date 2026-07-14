using MirrorWorld;

namespace MirrorWorldTests;

/// <summary>From-below pins for the adopted level-collision law (F129), the six-cosine
/// sibling of the seed world. Sources: docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md +
/// simulations/f129_level_collision_law.py (the committed gate carries the n &lt;= 210 census
/// and the named corner; these pins mirror its G2-G5 on the adopted range).</summary>
public class LevelCollisionTests
{
    static readonly World W = new();

    [Fact]
    public void TheLawHoldsOnTheAdoptedRange()
    {
        // G2-G4 shape: at every non-firing n <= 48 zero colliding pairs (mod-p distinctness
        // is PROOF of injectivity there); at every firing n at least one; sets equal.
        Assert.True(LevelCollision.LawHolds(5, 48));
    }

    [Theory]
    [InlineData(7)]
    [InlineData(11)]
    [InlineData(16)]
    [InlineData(25)]
    [InlineData(49)]
    public void NonFiringCombsAreInjective(int n)
    {
        Assert.False(LevelCollision.Fires(n));
        Assert.Equal(0, LevelCollision.CensusOf(n).CollidingPairs);
    }

    [Theory]
    [InlineData(9)]
    [InlineData(12)]
    [InlineData(20)]
    [InlineData(30)]
    public void FiringCombsCarryCollisions(int n)
    {
        Assert.True(LevelCollision.Fires(n));
        var census = LevelCollision.CensusOf(n);
        Assert.True(census.CollidingPairs > 0);
        Assert.NotNull(census.ExampleA);
    }

    [Fact]
    public void TheFiringPredicateHasItsThresholds()
    {
        // 3|n needs n >= 9 (n = 6 is silent); 10|n needs n >= 20 (n = 10 is silent).
        Assert.False(LevelCollision.Fires(6));
        Assert.False(LevelCollision.Fires(10));
        Assert.True(LevelCollision.Fires(9));
        Assert.True(LevelCollision.Fires(20));
        Assert.Equal(0, LevelCollision.CensusOf(6).CollidingPairs);
        Assert.Equal(0, LevelCollision.CensusOf(10).CollidingPairs);
    }

    [Fact]
    public void TheSubLawHolds_Overlap1Forces3DividesN()
    {
        Assert.True(LevelCollision.SubLawHolds(5, 48));
    }

    [Theory]
    [InlineData(20)]
    [InlineData(40)]
    public void ThePentagonDoorIsExclusivelyDisjoint(int n)
    {
        // 10|n with 3 not dividing n: collisions exist but every pair is fully disjoint.
        var census = LevelCollision.CensusOf(n);
        Assert.True(census.CollidingPairs > 0);
        Assert.Equal(0, census.Overlap1Pairs);
        Assert.Equal(census.CollidingPairs, census.DisjointPairs);
    }

    [Fact]
    public void TheMechanismAnchorsAreExact()
    {
        // n = 15: (8,12,14) ~ (9,11,13) at a nonzero shared level, four rotated R3 cycles;
        // n = 20: (1,7,9) ~ (3,5,10), the R5 conjugate pair + the zero mode.
        Assert.True(LevelCollision.AnchorsExact());
    }

    [Fact]
    public void TheCensusCountsMatchTheCommittedGateAtNamedCombs()
    {
        // pinned against the committed gate's exact layer (simulations/f129_level_collision_law.py
        // level_vec grouping, recomputed 2026-07-14: n=12 -> 25 exact pairs, n=20 -> 20).
        Assert.Equal(25, LevelCollision.CensusOf(12).CollidingPairs);
        Assert.Equal(20, LevelCollision.CensusOf(20).CollidingPairs);
    }

    [Fact]
    public void TheObjectOwnsItsTwoReadings()
    {
        var lc = new LevelCollision(W, 12);
        Assert.Equal(new[] { "levels", "collisions" }, lc.Own);
        Assert.Equal(12, lc.Census().N);
    }
}
