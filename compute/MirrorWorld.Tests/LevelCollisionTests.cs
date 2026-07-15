using MirrorWorld;

namespace MirrorWorldTests;

/// <summary>From-below pins for the adopted level-collision law (F129), the six-cosine
/// sibling of the seed world. Sources: docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md +
/// simulations/f129_level_collision_law.py (the committed gate carries the n &lt;= 210 census;
/// the named corner closed empty 2026-07-15; these pins mirror its G2-G5 on the adopted range).</summary>
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
    public void TheObjectOwnsItsThreeReadings()
    {
        var lc = new LevelCollision(W, 12);
        Assert.Equal(new[] { "levels", "collisions", "inventory" }, lc.Own);
        Assert.Equal(12, lc.Census().N);
    }

    // ---- the family inventory (adopted 2026-07-15): the thirteen derived closed forms tie
    // ---- to the independent GF(p) census exactly, total AND d-split, at every n.

    [Fact]
    public void TheInventoryMatchesTheCensusRowForRow()
    {
        // the from-below tie on the adopted range: the thirteen closed forms sum to the
        // census total, and the d-split (B/D/G/H share a mode, the rest disjoint) matches
        // too -- at FIRING and non-firing n alike (the inventory is 0 where the law is
        // injective; its formulas vanish at the silent doors n = 6, 10 by themselves).
        for (int n = 5; n <= 66; n++)
        {
            var census = LevelCollision.CensusOf(n);
            var inv = LevelCollision.InventoryOf(n);
            Assert.Equal(census.CollidingPairs, inv.Total);
            Assert.Equal(census.DisjointPairs, inv.Disjoint);
            Assert.Equal(census.Overlap1Pairs, inv.Overlap1);
        }
    }

    [Fact]
    public void TheCornerClosureSecondMechanismWalksAt70()
    {
        // n = 70: only the pentagon C and the corner-closure's second mechanism L
        // (zero + (R7:R5)) fire: 2*(70-10) = 120 disjoint + 20 disjoint, nothing shared.
        Assert.Equal(120, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.C, 70));
        Assert.Equal(20, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.L, 70));
        var census = LevelCollision.CensusOf(70);
        Assert.Equal(140, census.CollidingPairs);
        Assert.Equal(140, census.DisjointPairs);
        Assert.Equal(0, census.Overlap1Pairs);
        Assert.Equal((140L, 140L, 0L), LevelCollision.InventoryOf(70));
    }

    [Fact]
    public void TheOrderTwoTenDoorOpensAt105()
    {
        // n = 105, the first 105|n comb: seven families co-fire (doors 3, 15, 21, 105);
        // per-family closed forms and the census total 8858, tied exactly.
        Assert.Equal(5457, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.A, 105));
        Assert.Equal(1152, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.D, 105));
        Assert.Equal(612, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.E, 105));
        Assert.Equal(901, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.F, 105));
        Assert.Equal(576, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.H, 105));
        Assert.Equal(60, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.I, 105));
        Assert.Equal(100, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.M, 105));
        var census = LevelCollision.CensusOf(105);
        var inv = LevelCollision.InventoryOf(105);
        Assert.Equal(8858, inv.Total);
        Assert.Equal(census.CollidingPairs, inv.Total);
        Assert.Equal(census.DisjointPairs, inv.Disjoint);
        Assert.Equal(census.Overlap1Pairs, inv.Overlap1);
    }

    [Fact]
    public void TheOnsetsAreZerosOfTheInventory()
    {
        // the F129 thresholds are not extra conditions: the formulas vanish there themselves.
        Assert.Equal(0, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.A, 6));
        Assert.Equal(0, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.B, 6));
        Assert.Equal(0, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.C, 10));
        // ... and off-door every family is 0 (no door divides n = 25).
        foreach (LevelCollision.CollisionFamily f in Enum.GetValues<LevelCollision.CollisionFamily>())
            Assert.Equal(0, LevelCollision.FamilyCount(f, 25));
    }

    [Fact]
    public void TheParitySplitsOfEAndF()
    {
        // the counts proof section 5: the even-n deficits are single extra excluded labels
        // (E: c = 5u/2; F: the coarser special lattice), one formula per parity.
        Assert.Equal(12, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.E, 15));
        Assert.Equal(92, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.E, 30));
        Assert.Equal(212, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.E, 45));
        Assert.Equal(1, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.F, 15));
        Assert.Equal(25, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.F, 30));
        Assert.Equal(301, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.F, 45));
    }

    [Fact]
    public void TheMSplitAndTheImpossibleMiddleType()
    {
        // M sub-classifies 40 + 0 + 60 over the three CDK order-210 types; the middle type
        // (R7:R3,(R5:R3)) can NEVER fire (two distinct-type branches, one axis-fixed vertex).
        var (fanned, middle, fixedR5) = LevelCollision.MSplit;
        Assert.Equal(40, fanned);
        Assert.Equal(0, middle);
        Assert.Equal(60, fixedR5);
        Assert.Equal(LevelCollision.FamilyCount(LevelCollision.CollisionFamily.M, 105), fanned + middle + fixedR5);
    }

    [Fact]
    public void TheSubLawIsVisibleInTheDoors()
    {
        // every d = 2 family's door carries the factor 3 (overlap-1 forces 3|n, piece by piece).
        foreach (LevelCollision.CollisionFamily f in Enum.GetValues<LevelCollision.CollisionFamily>())
            if (LevelCollision.FamilySharesAMode(f))
                Assert.Equal(0, LevelCollision.FamilyDoor(f) % 3);
    }
}
