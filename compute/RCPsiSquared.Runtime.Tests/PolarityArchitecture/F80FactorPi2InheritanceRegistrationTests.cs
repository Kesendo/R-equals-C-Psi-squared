using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F80FactorPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop();

    [Fact]
    public void RegisterF80FactorPi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF80FactorPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F80FactorPi2Inheritance>());
    }

    [Fact]
    public void RegisterF80FactorPi2Inheritance_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF80FactorPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<F80FactorPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF80FactorPi2Inheritance_AncestorsContainBothPi2Axes()
    {
        // The composition: F80 inherits from BOTH Pi2-axes (number-anchor 2 from the
        // dyadic ladder, angle-anchor i from the Z₄ memory loop). Both must surface
        // in the registry's transitive ancestor closure.
        var registry = BuildBaseRegistry()
            .RegisterF80FactorPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F80FactorPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2I4MemoryLoopClaim), ancestors);
    }

    [Fact]
    public void RegisterF80FactorPi2Inheritance_CompositionAgreesAcrossTheRegistry()
    {
        // The "+2i" factor recomputed live through the registry must equal exactly
        // the parent ladder's Term(0) times the parent loop's PowerOfI(1).
        // Drift between the inheritance claim and either parent surfaces here.
        var registry = BuildBaseRegistry()
            .RegisterF80FactorPi2Inheritance()
            .Build();

        var f80 = registry.Get<F80FactorPi2Inheritance>();
        var ladder = registry.Get<Pi2DyadicLadderClaim>();
        var loop = registry.Get<Pi2I4MemoryLoopClaim>();

        Assert.Equal(ladder.Term(0), f80.TwoFactor, precision: 14);
        Assert.Equal(loop.PowerOfI(1), f80.IFactor);
        Assert.Equal(new Complex(0, 2), f80.PlusTwoIFactor);
        Assert.Equal(new Complex(0, -2), f80.MinusTwoIFactor);
    }

    [Fact]
    public void RegisterF80FactorPi2Inheritance_WithoutI4Loop_Throws()
    {
        // Without the Z₄ loop side, the inheritance edge cannot be drawn.
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: RegisterPi2I4MemoryLoop
                .RegisterF80FactorPi2Inheritance()
                .Build());
    }
}
