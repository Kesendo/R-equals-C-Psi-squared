using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class Pi2I4MemoryLoopRegistrationTests
{
    [Fact]
    public void RegisterPi2I4MemoryLoop_AddsClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop()
            .Build();

        Assert.True(registry.Contains<Pi2I4MemoryLoopClaim>());
    }

    [Fact]
    public void RegisterPi2I4MemoryLoop_TierIsTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<Pi2I4MemoryLoopClaim>().Tier);
    }

    [Fact]
    public void RegisterPi2I4MemoryLoop_AncestorsContainBothAxes()
    {
        // The dual-axis foundation: 90° angle-anchor (NinetyDegreeMirrorMemoryClaim)
        // AND the multiplicative inversion ladder (Pi2DyadicLadderClaim). Both must
        // surface as ancestors so the inheritance graph reads "Z₂ multiplicative ⊕ Z₄
        // rotational = complete mirror foundation."
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop()
            .Build();

        var ancestors = registry.AncestorsOf<Pi2I4MemoryLoopClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(NinetyDegreeMirrorMemoryClaim), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterPi2I4MemoryLoop_LiveMemoryClosureIsOne()
    {
        // The defining identity surfaces through the registry: i^4 = 1 + 0i exactly.
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop()
            .Build();

        var loop = registry.Get<Pi2I4MemoryLoopClaim>();
        Assert.Equal(new Complex(1, 0), loop.MemoryClosure());
    }

    [Fact]
    public void RegisterPi2I4MemoryLoop_WithoutDyadicLadder_Throws()
    {
        // Without the multiplicative-axis sibling, the dual-axis edge cannot be drawn.
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                // Missing: RegisterPi2DyadicLadder
                .RegisterPi2I4MemoryLoop()
                .Build());
    }
}
