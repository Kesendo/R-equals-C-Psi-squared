using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class Pi2Z4KleinDistinctionRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop();

    [Fact]
    public void RegisterPi2Z4KleinDistinction_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterPi2Z4KleinDistinction()
            .Build();

        Assert.True(registry.Contains<Pi2Z4KleinDistinctionClaim>());
    }

    [Fact]
    public void RegisterPi2Z4KleinDistinction_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterPi2Z4KleinDistinction()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<Pi2Z4KleinDistinctionClaim>().Tier);
    }

    [Fact]
    public void RegisterPi2Z4KleinDistinction_AncestorsContainBothSides()
    {
        // The distinction is between the two sides; the Object Manager makes this
        // queryable: AncestorsOf the distinction returns both Pi2I4MemoryLoop (Z₄)
        // and KleinFourCellClaim (Klein Z₂ × Z₂).
        var registry = BuildBaseRegistry()
            .RegisterPi2Z4KleinDistinction()
            .Build();

        var ancestors = registry.AncestorsOf<Pi2Z4KleinDistinctionClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2I4MemoryLoopClaim), ancestors);
        Assert.Contains(typeof(KleinFourCellClaim), ancestors);
    }

    [Fact]
    public void RegisterPi2Z4KleinDistinction_IsomorphismIsFalse()
    {
        // Tom's question answered through the registry: NOT isomorphic.
        var registry = BuildBaseRegistry()
            .RegisterPi2Z4KleinDistinction()
            .Build();

        var distinction = registry.Get<Pi2Z4KleinDistinctionClaim>();
        Assert.False(distinction.AreIsomorphic);
    }

    [Fact]
    public void RegisterPi2Z4KleinDistinction_SectorDimensionsMatchAtN3()
    {
        // The cardinality coincidence at N=3 (16-dim sectors on both sides) — the part
        // that tempted the conjecture. Verified through the registered claim.
        var registry = BuildBaseRegistry()
            .RegisterPi2Z4KleinDistinction()
            .Build();

        var distinction = registry.Get<Pi2Z4KleinDistinctionClaim>();
        Assert.Equal(16, distinction.SectorDimension(3));
    }

    [Fact]
    public void RegisterPi2Z4KleinDistinction_WithoutI4Loop_Throws()
    {
        // Without the Z₄ side, the distinction edge cannot be drawn.
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: RegisterPi2I4MemoryLoop
                .RegisterPi2Z4KleinDistinction()
                .Build());
    }
}
