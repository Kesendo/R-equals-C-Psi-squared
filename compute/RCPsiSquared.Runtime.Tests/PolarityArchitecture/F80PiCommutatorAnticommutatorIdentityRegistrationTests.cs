using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F80PiCommutatorAnticommutatorIdentityRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1PalindromeIdentity();

    [Fact]
    public void RegisterF80PiCommutatorAnticommutatorIdentity_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF80PiCommutatorAnticommutatorIdentity()
            .Build();

        Assert.True(registry.Contains<F80PiCommutatorAnticommutatorIdentity>());
    }

    [Fact]
    public void RegisterF80PiCommutatorAnticommutatorIdentity_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF80PiCommutatorAnticommutatorIdentity()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F80PiCommutatorAnticommutatorIdentity>().Tier);
    }

    [Fact]
    public void RegisterF80PiCommutatorAnticommutatorIdentity_AncestorsContainF1()
    {
        // The Π that conjugates the commutator is F1's Π; F1PalindromeIdentity is the
        // typed parent and must surface in the transitive ancestor closure.
        var registry = BuildBaseRegistry()
            .RegisterF80PiCommutatorAnticommutatorIdentity()
            .Build();

        var ancestors = registry.AncestorsOf<F80PiCommutatorAnticommutatorIdentity>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
    }

    [Fact]
    public void RegisterF80PiCommutatorAnticommutatorIdentity_SignTableHoldsAcrossRegistry()
    {
        // Cross-registry verification of the proven Step 5 sign table.
        var registry = BuildBaseRegistry()
            .RegisterF80PiCommutatorAnticommutatorIdentity()
            .Build();

        var f = registry.Get<F80PiCommutatorAnticommutatorIdentity>();
        Assert.Equal(1, f.SignFor('X', 'Y'));
        Assert.Equal(1, f.SignFor('Y', 'X'));
        Assert.Equal(-1, f.SignFor('X', 'Z'));
        Assert.Equal(-1, f.SignFor('Z', 'X'));
    }

    [Fact]
    public void RegisterF80PiCommutatorAnticommutatorIdentity_WithoutF1Parent_Throws()
    {
        // The F1 parent edge is required; without F1PalindromeIdentity registered the
        // inheritance edge cannot be drawn.
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF80PiCommutatorAnticommutatorIdentity()
                .Build());
    }
}
